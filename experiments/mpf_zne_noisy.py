"""MPF + ZNE under realistic noise: does combining them beat the previous report's approach?

For each time t on a short-time grid, on one fixed zero-SWAP layout of a fake IBM backend:
  * Trotter circuits with k steps (k in the MPF set, max k = k_max) are folded at scales 1,3,5
  * every circuit is run with Aer + the backend's noise model
Methods compared against exact ED (per observable):
  raw          Trotter k_max, no mitigation
  linZNE       Trotter k_max + linear ZNE              (the previous report's method)
  expZNE       Trotter k_max + exponential ZNE
  richZNE      Trotter k_max + Richardson ZNE
  MPF+expZNE   exponential ZNE on each k, then MPF combination   (ours)
  MPF+richZNE  Richardson ZNE on each k, then MPF combination
plus the noiseless references ideal Trotter (k_max) and ideal MPF.
Error metric: the report's, 100 * max_t |dO| / max_t |O_ED| over the grid, as mean +- std over
--reps independent simulation seeds (one seed is too noisy to rank methods); also |dO| per t.

  python experiments/mpf_zne_noisy.py                       # N=6 chain, h/J = 0.5 and 2
  python experiments/mpf_zne_noisy.py --h 1 --ks 4,6,8 --l1 1.5 --shots 8000
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qfest import hardware as hw
from qfest.ed import quench_state
from qfest.extrapolation import exponential, richardson, linear
from qfest.metrics import max_report_deviation
from qfest.mpf import approx_coefficients, noise_amplification, richardson_coefficients
from qfest.noise import simulator
from qfest.results import RESULTS_DIR, save
from qfest.shots import counts_expectations, state_expectations
from qfest.tfim import tfim_circuit, trotter_state

OBS = ("Mzz", "Zavg")
# Recommended MPF sets from experiments/mpf_short_times.py (shot noise only):
#   h/J <= 1: [4,6,8] with ||x||_1 <= 1.5 ; h/J = 2: exact [6,8]
DEFAULT_MPF = {0.5: ([4, 6, 8], 1.5), 1.0: ([4, 6, 8], 1.5), 2.0: ([6, 8], None)}
BLUE, ORANGE, AQUA, YELLOW, INK, INK2, GRID, SURFACE = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb")


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--ring", action="store_true", help="periodic ring (needs n=12 on heavy-hex; slow)")
    p.add_argument("--h", type=float, nargs="+", default=[0.5, 2.0])
    p.add_argument("--tmax", type=float, default=2.0)
    p.add_argument("--nt", type=int, default=8)
    p.add_argument("--ks", default=None, help="MPF step counts, e.g. 4,6,8 (default: per-h recommendation)")
    p.add_argument("--l1", type=float, default=None, help="||x||_1 bound for MPF coefficients")
    p.add_argument("--scales", type=int, nargs="+", default=[1, 3, 5])
    p.add_argument("--shots", type=int, default=8000, help="shots per circuit")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--reps", type=int, default=5, help="independent shot/noise seeds (seed, seed+1, ...)")
    p.add_argument("--fake", default="marrakesh")
    p.add_argument("--tag", default="", help="suffix for output file names")
    p.add_argument("--no-plot", action="store_true")
    return p.parse_args()


def mpf_setup(a, h):
    if a.ks:
        ks, l1 = [int(k) for k in a.ks.split(",")], a.l1
    else:
        ks, l1 = DEFAULT_MPF.get(h, ([4, 6, 8], 1.5))
    x = approx_coefficients(ks, l1) if l1 else richardson_coefficients(ks)
    return ks, x, l1


def run_h(a, h, backend, layout, sim):
    n, periodic = a.n, a.ring
    ks, x, l1 = mpf_setup(a, h)
    kmax = max(ks)
    times = np.linspace(a.tmax / a.nt, a.tmax, a.nt)

    jobs, index = [], {}
    for i, t in enumerate(times):
        for k in ks:
            isa = hw.to_isa(tfim_circuit(n, 1.0, h, t / k, k, periodic=periodic), backend, layout, 1, a.seed)
            circs, factors = hw.zne_circuits(isa, n, backend, a.scales)
            index[(i, k)] = (len(jobs), factors)
            jobs.extend(circs)
    reps = [analyse(a, h, n, periodic, ks, x, kmax, times, index,
                    hw.run_simulator(jobs, sim, a.shots, a.seed + r)) for r in range(a.reps)]
    first = reps[0]
    names = list(first["summary"][OBS[0]].keys())
    summary = {o: {m: [rep["summary"][o][m] for rep in reps] for m in names} for o in OBS}
    stats = {o: {m: (float(np.mean(v)), float(np.std(v))) for m, v in summary[o].items()} for o in OBS}
    return {"h": h, "times": times, "ks": ks, "coefficients": x.tolist(), "l1_bound": l1,
            "noise_amplification": noise_amplification(x), "exact": first["exact"],
            "ideal_trotter": first["ideal"], "ideal_mpf": first["ideal_mpf"], "methods": first["methods"],
            "report_dev_per_rep": summary, "report_dev_mean_std": stats,
            "circuits_per_t": len(ks) * len(a.scales), "shots_per_circuit": a.shots, "reps": a.reps}


def analyse(a, h, n, periodic, ks, x, kmax, times, index, counts):
    exact = {o: [] for o in OBS}
    ideal = {o: [] for o in OBS}
    ideal_mpf = {o: [] for o in OBS}
    methods = {m: {o: [] for o in OBS} for m in ("raw", "linZNE", "expZNE", "richZNE", "MPF+expZNE", "MPF+richZNE")}
    for i, t in enumerate(times):
        ex = state_expectations(quench_state(n, t, 1.0, h, periodic), n, periodic)
        ideal_k = {k: state_expectations(trotter_state(n, 1.0, h, t, k, periodic=periodic), n, periodic) for k in ks}
        raw = {}
        for k in ks:
            start, factors = index[(i, k)]
            raw[k] = [counts_expectations(counts[start + j], n, periodic) for j in range(len(a.scales))]
        for o in OBS:
            exact[o].append(ex[o])
            ideal[o].append(ideal_k[kmax][o])
            ideal_mpf[o].append(float(sum(c * ideal_k[k][o] for c, k in zip(x, ks))))
            f = index[(i, kmax)][1]
            vals = [r[o] for r in raw[kmax]]
            methods["raw"][o].append(vals[0])
            methods["linZNE"][o].append(linear(f, vals))
            methods["expZNE"][o].append(exponential(f, vals))
            methods["richZNE"][o].append(richardson(f, vals))
            for name, fit in (("MPF+expZNE", exponential), ("MPF+richZNE", richardson)):
                per_k = [fit(index[(i, k)][1], [r[o] for r in raw[k]]) for k in ks]
                methods[name][o].append(float(sum(c * v for c, v in zip(x, per_k))))

    summary = {o: {m: max_report_deviation(v[o], exact[o]) for m, v in
                   {**methods, "ideal_trotter": ideal, "ideal_MPF": ideal_mpf}.items()} for o in OBS}
    return {"exact": exact, "ideal": ideal, "ideal_mpf": ideal_mpf, "methods": methods, "summary": summary}


def plot(results, out, obs="Mzz"):
    fig, axes = plt.subplots(2, len(results), figsize=(5.2 * len(results), 7.0), dpi=150, squeeze=False)
    series = (("raw", "raw (no mitigation)", INK2, "o"), ("linZNE", "Trotter + linear ZNE (old report)", ORANGE, "s"),
              ("expZNE", "Trotter + exponential ZNE", BLUE, "^"), ("MPF+expZNE", "MPF + exponential ZNE (ours)", AQUA, "D"))
    for col, r in enumerate(results):
        t = r["times"]
        top, bot = axes[0, col], axes[1, col]
        top.plot(t, r["exact"][obs], color=INK, linewidth=2, label="exact (ED)")
        top.plot(t, r["ideal_trotter"][obs], color=INK2, linewidth=1.2, linestyle="--", label="ideal Trotter (noiseless)")
        for key, label, color, marker in series:
            top.plot(t, r["methods"][key][obs], color=color, linewidth=2 if key != "raw" else 1.2,
                     marker=marker, markersize=7, label=label)
            err = np.abs(np.array(r["methods"][key][obs]) - np.array(r["exact"][obs]))
            bot.plot(t, err, color=color, linewidth=2 if key != "raw" else 1.2, marker=marker, markersize=7, label=label)
            mean, std = r["report_dev_mean_std"][obs][key]
            bot.annotate(f"{mean:.1f}±{std:.1f}%", (t[-1], err[-1]), textcoords="offset points",
                         xytext=(6, 0), va="center", fontsize=8, color=INK2)
        for ax, title, ylabel in ((top, f"{obs} vs time, h/J = {r['h']:g}", obs),
                                  (bot, f"|error| vs exact, seed 1 (labels: report % mean±std over seeds)", f"|Δ{obs}|")):
            ax.set_title(title, color=INK, fontsize=10, loc="left")
            ax.set_xlabel("time t (1/J)", color=INK2)
            ax.set_ylabel(ylabel, color=INK2)
            ax.grid(True, color=GRID, linewidth=0.8)
            ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.tick_params(colors=INK2)
        top_err = max(np.max(np.abs(np.array(r["methods"][k][obs]) - np.array(r["exact"][obs]))) for k, *_ in series)
        bot.set_ylim(-0.06 * top_err, 1.08 * top_err)
        bot.set_xlim(right=t[-1] * 1.22)
    axes[0, 0].legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle(f"MPF + ZNE under FakeMarrakesh noise (N=6 chain, zero-SWAP layout, k_max = 8, {results[0]['shots_per_circuit']} shots/circuit)",
                 color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)


def replot(paths, out):
    """Redraw the figure from saved results JSON files (no simulation)."""
    import json
    plot([json.loads(open(p).read())["extra"] for p in paths], out)


def main():
    a = parse()
    backend = hw.get_backend(fake=a.fake)
    layout = hw.find_layout(backend, a.n, a.ring)
    sim = simulator(backend, "backend")
    results = []
    for h in a.h:
        r = run_h(a, h, backend, layout, sim)
        results.append(r)
        print(f"\nh/J={h}  MPF ks={r['ks']} l1<={r['l1_bound']}  ||x||1={r['noise_amplification']:.2f}  "
              f"({r['circuits_per_t']} circuits x {a.shots} shots per time point)")
        print(f"{'method':<14}" + "".join(f"{o:>16}" for o in OBS) + f"   (report %, mean ± std over {a.reps} seeds)")
        for m in ("ideal_trotter", "ideal_MPF", "raw", "linZNE", "richZNE", "expZNE", "MPF+richZNE", "MPF+expZNE"):
            print(f"{m:<14}" + "".join(f"{r['report_dev_mean_std'][o][m][0]:>9.2f} ±{r['report_dev_mean_std'][o][m][1]:<5.2f}"
                                       for o in OBS))
        save({"experiment": "mpf_zne_noisy", "backend": f"aer_backend_fake_{a.fake}",
              "params": {**vars(a), "h": h, "layout": layout}, "method": "mpf_zne",
              "t": r["times"], "observables": {"exact": r["exact"], "ideal_trotter": r["ideal_trotter"]},
              "extra": r}, f"mpf_zne_noisy_n{a.n}_{'ring' if a.ring else 'chain'}_h{h}{a.tag}")
    if a.no_plot:
        return
    figs = RESULTS_DIR / "figs"
    figs.mkdir(parents=True, exist_ok=True)
    out = figs / f"mpf_zne_noisy_n{a.n}_{'ring' if a.ring else 'chain'}{a.tag}.png"
    plot(results, out)
    print(f"\nfigure: {out}")


if __name__ == "__main__":
    main()
