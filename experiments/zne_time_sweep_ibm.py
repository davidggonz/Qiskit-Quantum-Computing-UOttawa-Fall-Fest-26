"""ZNE vs time on a real IBM QPU (hardware twin of experiments/zne_time_sweep.py).

Three steps, so the QPU job and the analysis are decoupled:

  plan     transpile everything, print 2q-gate counts and depth per circuit (no QPU time)
  submit   ONE EstimatorV2 job with every (time, noise factor) circuit; saves the job id
  analyze  fetch the job by id, extrapolate to zero noise, compare with ED, plot

Physics: N = 12 periodic ring (0-SWAP heavy-hex embedding), quench from |0...0>,
2nd-order Trotter (dt = 0.1), observables Mz (RMS), Mx, Mzz. ZNE by local folding of
every 2q gate (G -> G (G† G)^k, lambda = 1, 3, 5) on the transpiled circuit, then the
team's extrapolators (linear, richardson, exponential) plus exp_fixed: exponential with
its asymptote pinned to the fully depolarized value (Mz -> sqrt(1/N), Mx, Mzz -> 0). Runtime adds TREX readout
mitigation, Pauli twirling and dynamical decoupling; Runtime's own ZNE stays OFF so the
folding is ours and directly comparable with the simulation.

Examples:
  python experiments/zne_time_sweep_ibm.py plan
  python experiments/zne_time_sweep_ibm.py submit --backend ibm_fez
  python experiments/zne_time_sweep_ibm.py analyze --job-id <id>
  python experiments/zne_time_sweep_ibm.py analyze --from-json results/zne_time_sweep_ibm_h1_ibm_quebec.json
  python experiments/zne_time_sweep_ibm.py submit --fake      # local test, no token, no QPU
"""
import argparse
import json
import pathlib
import sys
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import hardware as hw  # noqa: E402
from qfest.ed import quench  # noqa: E402
from qfest.extrapolation import EXTRAPOLATORS  # noqa: E402
from qfest.results import save  # noqa: E402

N, J, DT = 12, 1.0, 0.1
LAMS = [1, 3, 5]
METHODS = ["linear", "richardson", "exponential", "exp_fixed"]
# Fully depolarized limits: (Σz)²/N² -> 1/N so Mz -> sqrt(1/N); <X>, <ZZ> -> 0.
# exp_fixed pins the exponential's asymptote there: 2 parameters from 3 noise
# factors, so a goodness-of-fit (chi2, 1 dof) exists, unlike the free 3-parameter fit.
ASYMPTOTE = {"Mz": float(np.sqrt(1.0 / N)), "Mx": 0.0, "Mzz": 0.0}
LABEL = {"raw": "raw (λ=1)", "linear": "linear", "richardson": "richardson",
         "exponential": "exponential (free asymptote)", "exp_fixed": "exponential (fixed asymptote)"}
KEYS = ["Mz", "Mx", "Mzz"]
EST_KEYS = ["Mz2", "Mx", "Mzz"]
JOBS_DIR = ROOT / "results" / "ibm_jobs"
def exact_at(times, h):
    """ED observables at arbitrary (not necessarily uniform) times."""
    pts = [quench(N, [t], J, h) for t in times]
    return {k: np.array([p[k][0] for p in pts]) for k in KEYS}


COL = {"raw": "#52514e", "linear": "#eb6834", "richardson": "#1baf7a", "exponential": "#2a78d6",
       "exp_fixed": "#4a3aa7"}


def extrapolate(method, k, vals):
    if method == "exp_fixed":
        return float(EXTRAPOLATORS["exponential"](LAMS, vals, asymptote=ASYMPTOTE[k]))
    return float(EXTRAPOLATORS[method](LAMS, vals))


def fixed_exp_chi2(k, vals, errs):
    """chi2 (1 dof) of a*exp(-b*lam) + asymptote through the 3 noise factors."""
    from scipy.optimize import curve_fit
    lam = np.asarray(LAMS, float)
    f = lambda x, a, b: a * np.exp(-b * x) + ASYMPTOTE[k]
    try:
        p, _ = curve_fit(f, lam, vals, p0=[vals[0] - ASYMPTOTE[k], 0.3], maxfev=10000)
    except (RuntimeError, ValueError):
        return float("nan")
    return float(np.sum(((vals - f(lam, *p)) / np.maximum(errs, 1e-9)) ** 2))


def get_backend(args):
    if args.fake:
        print("[backend] FakeMarrakesh in local mode (no QPU, no token)")
        return hw.get_fake_backend("marrakesh")
    service = hw.get_service()
    if service is None:
        sys.exit("No QISKIT_IBM_TOKEN: copy .env.example to .env and paste your token, or use --fake")
    return service.backend(args.backend) if args.backend else hw.select_backend(service)


def build_pubs(backend, h, times):
    """One PUB per (t, lambda): folded ISA circuit + the 3 observables mapped to its layout."""
    from qiskit import transpile
    layout = hw.select_best_qubits(backend, N, periodic=True)
    obs = hw.build_tfim_observables(N, J, h, periodic=True)
    pubs, index = [], []
    for t in times:
        steps = int(round(t / DT))
        base = transpile(hw.build_tfim_circuit(N, J, h, DT, steps, add_measure=False), backend,
                         optimization_level=1, initial_layout=layout, seed_transpiler=42)
        isa_obs = [obs[k].apply_layout(base.layout) for k in EST_KEYS]
        for lam in LAMS:
            qc = hw.fold_two_qubit_gates(base, lam)
            pubs.append((qc, isa_obs))
            ops = qc.count_ops()
            index.append({"t": t, "steps": steps, "lambda": lam,
                          "two_qubit_gates": hw.two_qubit_count(ops), "depth": qc.depth()})
    return pubs, index, layout


def cmd_plan(args):
    backend = get_backend(args)
    hw.backend_properties_dict(backend)
    pubs, index, layout = build_pubs(backend, args.h, args.times)
    print(f"\nlayout (12-ring): {layout}")
    print(f"{'t':>5} {'lam':>4} {'2q gates':>9} {'depth':>6}")
    for r in index:
        print(f"{r['t']:>5} {r['lambda']:>4} {r['two_qubit_gates']:>9} {r['depth']:>6}")
    print(f"\n{len(pubs)} circuits x {args.shots} shots, 3 observables each, in ONE job.")
    print("Rule of thumb: above ~1500 two-qubit gates the raw signal is close to noise.")


def cmd_submit(args):
    from qiskit_ibm_runtime import EstimatorV2
    backend = get_backend(args)
    pubs, index, layout = build_pubs(backend, args.h, args.times)
    opts = {
        "default_shots": args.shots,
        "resilience_level": 1,                       # TREX readout mitigation
        "twirling": {"enable_gates": True, "num_randomizations": "auto"},
        "dynamical_decoupling": {"enable": True, "sequence_type": "XpXm"},
    }
    if args.fake:  # local mode ignores runtime-only error suppression
        opts = {"default_shots": args.shots}
    est = EstimatorV2(mode=backend, options=opts)
    job = est.run(pubs)
    job_id = job.job_id()
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    meta = {"job_id": job_id, "backend": backend.name, "fake": args.fake, "h": args.h,
            "times": args.times, "lambdas": LAMS, "shots": args.shots, "layout": layout,
            "index": index, "options": opts,
            "submitted": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    path = JOBS_DIR / f"zne_sweep_h{args.h:g}_{job_id}.json"
    path.write_text(json.dumps(meta, indent=2))
    print(f"[submit] job id: {job_id}  (metadata saved to {path})")
    if args.fake:  # local mode finishes immediately: keep the result for `analyze --fake`
        res = job.result()
        cache = {"evs": [list(map(float, r.data.evs)) for r in res],
                 "stds": [list(map(float, r.data.stds)) for r in res]}
        (JOBS_DIR / f"{job_id}_result.json").write_text(json.dumps(cache))
    else:
        print("Check progress at https://quantum.cloud.ibm.com/workloads, then run:")
        print(f"  python experiments/zne_time_sweep_ibm.py analyze --job-id {job_id}")


def load_result(job_id, meta):
    cache = JOBS_DIR / f"{job_id}_result.json"
    if cache.exists():
        d = json.loads(cache.read_text())
        return np.array(d["evs"]), np.array(d["stds"])
    service = hw.get_service()
    job = service.job(job_id)
    print(f"[analyze] status: {job.status()}")
    res = job.result()
    evs = np.array([list(map(float, r.data.evs)) for r in res])
    stds = np.array([list(map(float, r.data.stds)) for r in res])
    cache.write_text(json.dumps({"evs": evs.tolist(), "stds": stds.tolist()}))
    return evs, stds


def load_from_results_json(path):
    """Rebuild (meta, evs, stds) from a results/zne_time_sweep_ibm_*.json written by analyze,
    so the analysis can be redone without the job cache or QPU access."""
    d = json.loads(pathlib.Path(path).read_text())
    ex = d["extra"]
    is_ibm = d["backend"].startswith("ibm_")
    meta = {"backend": d["backend"][len("ibm_"):] if is_ibm else "fake_local", "fake": not is_ibm,
            "h": d["params"]["h"], "times": d["t"], "shots": d["params"]["shots"],
            "layout": d["params"]["layout"], "index": ex["index"], "options": ex["options"],
            "job_id": ex["job_id"]}
    return meta, np.array(ex["raw_evs"]), np.array(ex["raw_stds"])


def cmd_analyze(args):
    if args.from_json:
        meta, evs, stds = load_from_results_json(args.from_json)
        args.job_id = meta["job_id"]
    else:
        metas = sorted(JOBS_DIR.glob(f"zne_sweep_*_{args.job_id}.json"))
        if not metas:
            sys.exit(f"no metadata for job {args.job_id} in {JOBS_DIR} (or pass --from-json)")
        meta = json.loads(metas[0].read_text())
        evs, stds = load_result(args.job_id, meta)
    times, h = meta["times"], meta["h"]
    exact = exact_at(times, h)
    rng = np.random.default_rng(0)

    def to_obs(ev):  # columns Mz2, Mx, Mzz -> Mz, Mx, Mzz
        return np.stack([np.sqrt(np.clip(ev[..., 0], 0, None)), ev[..., 1], ev[..., 2]], axis=-1)

    # Monte Carlo over the reported standard errors -> error bars on every estimator
    n_mc = 300
    agg = {m: {k: {"mean": [], "std": [], "abs_err": []} for k in KEYS} for m in ["raw"] + METHODS}
    chi2 = {k: [] for k in KEYS}
    for i, t in enumerate(times):
        rows = [j for j, r in enumerate(meta["index"]) if r["t"] == t]
        samples = to_obs(evs[rows][None] + stds[rows][None] * rng.standard_normal((n_mc, len(rows), 3)))
        central = to_obs(evs[rows])
        for c, k in enumerate(KEYS):
            ests = {"raw": (central[0, c], samples[:, 0, c])}
            for m in METHODS:
                ests[m] = (extrapolate(m, k, central[:, c]),
                           np.array([extrapolate(m, k, s[:, c]) for s in samples]))
            chi2[k].append(fixed_exp_chi2(k, central[:, c], samples[:, :, c].std(axis=0)))
            for m, (val, mc) in ests.items():
                agg[m][k]["mean"].append(float(val))
                agg[m][k]["std"].append(float(np.std(mc)))
                agg[m][k]["abs_err"].append(float(abs(val - exact[k][i])))

    tag = f"h{h:g}_{meta['backend']}"
    save({
        "experiment": f"zne_time_sweep_ibm_{tag}",
        "backend": f"ibm_{meta['backend']}" if not meta["fake"] else "fake_local",
        "params": {"n": N, "J": J, "h": h, "dt": DT, "order": 2, "periodic": True,
                   "layout": meta["layout"], "noise_factors": LAMS, "shots": meta["shots"]},
        "method": "zne_" + "|".join(METHODS),
        "t": times,
        "observables": {"exact": {k: exact[k].tolist() for k in KEYS}, "estimators": agg},
        "extra": {"job_id": args.job_id, "index": meta["index"], "options": meta["options"],
                  "asymptotes_exp_fixed": ASYMPTOTE, "chi2_exp_fixed_1dof": chi2,
                  "raw_evs": evs.tolist(), "raw_stds": stds.tolist()},
    }, f"zne_time_sweep_ibm_{tag}")

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    for c, k in enumerate(KEYS):
        ax = axes[0, c]
        ax.plot(times, exact[k], color="#0b0b0b", lw=2, label="exact (ED)")
        for m in ["raw"] + METHODS:
            ax.errorbar(times, agg[m][k]["mean"], yerr=agg[m][k]["std"], color=COL[m], marker="o",
                        ms=5, lw=1.5, capsize=3, label=LABEL[m])
        ax.set_title(k)
        ax = axes[1, c]
        for m in ["raw"] + METHODS:
            ax.semilogy(times, np.maximum(agg[m][k]["abs_err"], 1e-5), color=COL[m], marker="o", ms=5, lw=1.5)
        ax.set_xlabel("time t (1/J)")
        for a in axes[:, c]:
            a.axvline(1.5, color="#a3a29d", ls=":", lw=1)
            a.grid(color="#e4e3df")
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
    axes[0, 0].set_ylabel("value")
    axes[1, 0].set_ylabel("|estimate − ED|")
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(f"ZNE vs time on {meta['backend']} (job {args.job_id[:12]}…), N = {N} ring, h/J = {h:g}, "
                 f"dt = {DT}, {meta['shots']} shots, TREX + twirling + DD, λ = 1, 3, 5\n"
                 f"exponential (fixed asymptote): Mz → √(1/N), Mx, Mzz → 0")
    fig.tight_layout()
    out = ROOT / "results" / "figs"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"zne_time_sweep_ibm_{tag}.png"
    fig.savefig(p, dpi=150)
    print(f"[analyze] figure: {p}")
    print(f"\n{'method':>12} {'mean|err|':>10}  per observable: |err| (sigma) at each t")
    for m in ["raw"] + METHODS:
        mean_err = np.mean([e for k in KEYS for e in agg[m][k]["abs_err"]])
        cells = {k: [f"{e:.3f}({e / max(s, 1e-12):.1f}σ)" for e, s in zip(agg[m][k]["abs_err"], agg[m][k]["std"])]
                 for k in KEYS}
        print(f"{m:>12} {mean_err:>10.4f}  {cells}")
    print("chi2 (1 dof) of the fixed-asymptote fit:", {k: [round(x, 1) for x in v] for k, v in chi2.items()})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["plan", "submit", "analyze"])
    ap.add_argument("--backend", default="", help="e.g. ibm_fez; default: preferred/least busy")
    ap.add_argument("--h", type=float, default=1.0)
    ap.add_argument("--times", default="0.5,1.0,1.5,2.0",
                    type=lambda s: [float(x) for x in s.split(",")])
    ap.add_argument("--shots", type=int, default=4000)
    ap.add_argument("--job-id", default="")
    ap.add_argument("--fake", action="store_true", help="FakeMarrakesh local mode (testing)")
    ap.add_argument("--from-json", default="",
                    help="analyze: re-analyze a saved results/zne_time_sweep_ibm_*.json (no QPU/cache needed)")
    args = ap.parse_args()
    {"plan": cmd_plan, "submit": cmd_submit, "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    main()
