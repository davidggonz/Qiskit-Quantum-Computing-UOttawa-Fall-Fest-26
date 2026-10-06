"""ZNE vs time on a realistic noise model: does exponential ZNE survive past t ~ 1.5,
where the previous report's linear ZNE started to overcorrect?

Setup: N = 12 periodic ring (0-SWAP heavy-hex embedding), quench from |0...0>,
2nd-order Trotter with dt = 0.1, noise = NoiseModel.from_backend(FakeMarrakesh)
on the actual ring qubits (per-qubit/per-coupler errors + thermal relaxation).
Readout error is left out on purpose: on hardware it is handled by TREX, and
gate folding does not amplify it, so ZNE cannot remove it.

Observables are estimated from shots (Z-basis circuit for Mz, Mzz; X-basis circuit
for Mx). Folding: every CZ -> CZ (CZ CZ)^k after transpiling, lambda = 1, 3, 5.
Error bars: std over `seeds` independent groups of `shots` shots (shot noise).

Run:  python experiments/zne_time_sweep.py [--h 1.0] [--seeds 3] [--shots 4000]
Writes results/zne_time_sweep_h<h>.json (team schema) and results/figs/zne_time_sweep_h<h>.png
"""
import argparse
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
METHODS = ["linear", "richardson", "exponential"]
KEYS = ["Mz", "Mx", "Mzz"]
def exact_at(times, h):
    """ED observables at arbitrary (not necessarily uniform) times."""
    pts = [quench(N, [t], J, h) for t in times]
    return {k: np.array([p[k][0] for p in pts]) for k in KEYS}


COL = {"raw": "#52514e", "linear": "#eb6834", "richardson": "#1baf7a", "exponential": "#2a78d6"}


def measure_point(sim, be, layout, h, steps, shots, seeds):
    """Raw noisy observables at each folding factor, one dict per seed.

    The density-matrix simulation is deterministic, so each circuit runs once with
    shots*seeds shots and the per-shot memory is split into `seeds` independent groups
    (error bars = shot noise at `shots` shots per circuit).
    """
    from qiskit import transpile
    vals = [{k: [] for k in KEYS} for _ in range(seeds)]
    circs = {b: transpile(hw.build_tfim_circuit(N, J, h, DT, steps, basis=b), be,
                          optimization_level=1, initial_layout=layout, seed_transpiler=42)
             for b in ("z", "x")}
    for lam in LAMS:
        groups = {}
        for b, qc in circs.items():
            mem = sim.run(hw.fold_two_qubit_gates(qc, lam), shots=shots * seeds, memory=True,
                          seed_simulator=17 + lam).result().get_memory()
            groups[b] = [mem[i * shots:(i + 1) * shots] for i in range(seeds)]
        for i in range(seeds):
            cz, cx = ({}, {})
            for bits in groups["z"][i]:
                cz[bits] = cz.get(bits, 0) + 1
            for bits in groups["x"][i]:
                cx[bits] = cx.get(bits, 0) + 1
            z = hw.observables_from_counts(cz, N)
            vals[i]["Mz"].append(z["Mz"])
            vals[i]["Mzz"].append(z["Mzz"])
            vals[i]["Mx"].append(hw.mx_from_counts(cx, N))
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=1.0)
    ap.add_argument("--times", default="0.5,1.0,1.5,2.0,2.5,3.0")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--shots", type=int, default=4000)
    args = ap.parse_args()
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel

    be = hw.get_fake_backend("marrakesh")
    layout = hw.select_best_qubits(be, N, periodic=True)
    sim = AerSimulator(method="density_matrix",
                       noise_model=NoiseModel.from_backend(be, readout_error=False))
    times = [float(t) for t in args.times.split(",")]
    exact = exact_at(times, args.h)
    rows = []
    t0 = time.time()
    for t in times:
        steps = int(round(t / DT))
        per_seed = []
        for raw in measure_point(sim, be, layout, args.h, steps, args.shots, args.seeds):
            est = {"raw": {k: raw[k][0] for k in KEYS}}
            for m in METHODS:
                est[m] = {k: float(EXTRAPOLATORS[m](LAMS, raw[k])) for k in KEYS}
            per_seed.append({"raw_by_lambda": raw, "estimates": est})
        rows.append({"t": t, "steps": steps, "seeds": per_seed})
        print(f"[t={t}] done ({time.time() - t0:.0f}s)", flush=True)

    # aggregate: mean and std over seeds of each estimator, error vs ED
    agg = {m: {k: {"mean": [], "std": [], "abs_err": []} for k in KEYS} for m in ["raw"] + METHODS}
    for i, r in enumerate(rows):
        for m in agg:
            for k in KEYS:
                v = np.array([s["estimates"][m][k] for s in r["seeds"]])
                agg[m][k]["mean"].append(float(v.mean()))
                agg[m][k]["std"].append(float(v.std(ddof=1)) if len(v) > 1 else 0.0)
                agg[m][k]["abs_err"].append(float(abs(v.mean() - exact[k][i])))

    tag = f"h{args.h:g}"
    save({
        "experiment": f"zne_time_sweep_{tag}",
        "backend": "aer_noisy (NoiseModel.from_backend(FakeMarrakesh), no readout error)",
        "params": {"n": N, "J": J, "h": args.h, "dt": DT, "order": 2, "periodic": True,
                   "layout": layout, "noise_factors": LAMS, "shots": args.shots, "seeds": args.seeds},
        "method": "zne_" + "|".join(METHODS),
        "t": times,
        "observables": {"exact": {k: exact[k].tolist() for k in KEYS}, "estimators": agg},
        "extra": {"per_point": rows, "wall_s": round(time.time() - t0, 1)},
    }, f"zne_time_sweep_{tag}")

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2]})
    for c, k in enumerate(KEYS):
        ax = axes[0, c]
        ax.plot(times, exact[k], color="#0b0b0b", lw=2, label="exact (ED)")
        for m in ["raw"] + METHODS:
            ax.errorbar(times, agg[m][k]["mean"], yerr=agg[m][k]["std"], color=COL[m],
                        marker="o", ms=5, lw=1.5, capsize=3, label=m if m != "raw" else "raw (λ=1)")
        ax.axvline(1.5, color="#a3a29d", ls=":", lw=1)
        ax.set_title(k)
        ax.grid(color="#e4e3df")
        ax = axes[1, c]
        for m in ["raw"] + METHODS:
            ax.semilogy(times, np.maximum(agg[m][k]["abs_err"], 1e-5), color=COL[m], marker="o", ms=5, lw=1.5)
        ax.axvline(1.5, color="#a3a29d", ls=":", lw=1)
        ax.set_xlabel("time t (1/J)")
        ax.grid(color="#e4e3df")
    axes[0, 0].set_ylabel("value")
    axes[1, 0].set_ylabel("|estimate − ED|")
    axes[0, 0].legend(fontsize=8, frameon=False)
    axes[0, 0].text(1.52, axes[0, 0].get_ylim()[0], " previous linear ZNE\n broke here", fontsize=7, color="#52514e",
                    va="bottom")
    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(f"ZNE vs time, N = {N} ring, h/J = {args.h:g}, dt = {DT}, FakeMarrakesh noise "
                 f"(no readout), {args.shots} shots × {args.seeds} seeds, λ = 1, 3, 5")
    fig.tight_layout()
    out = ROOT / "results" / "figs"
    out.mkdir(parents=True, exist_ok=True)
    p = out / f"zne_time_sweep_{tag}.png"
    fig.savefig(p, dpi=150)
    print(p)
    for m in ["raw"] + METHODS:
        print(m, {k: [round(e, 3) for e in agg[m][k]["abs_err"]] for k in KEYS})


if __name__ == "__main__":
    main()
