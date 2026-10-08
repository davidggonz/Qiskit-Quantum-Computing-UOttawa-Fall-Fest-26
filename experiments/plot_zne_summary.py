"""Summary figures for the IBM hardware ZNE sweeps (results section + slides).

Reads the analyzed hardware runs in reported_results/ and writes to reported_results/figs/:
  zne_summary.png           3 panels: (a) mean error per method and run,
                            (b) h/J = 1 error vs time, (c) h/J = 2 error vs time
  zne_summary_methods.png   panel (a) alone, slide-sized
  zne_summary_vs_time.png   panels (b, c), slide-sized

"Error" is |estimate - exact diagonalization|, averaged over Mz, Mx, Mzz; error bars
combine the per-observable shot-noise standard deviations.

Run:  python experiments/plot_zne_summary.py
"""
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES, FIGS = ROOT / "reported_results", ROOT / "reported_results" / "figs"
KEYS = ["Mz", "Mx", "Mzz"]

RUNS = [  # (file, label, marker)
    ("zne_time_sweep_ibm_h1_ibm_quebec.json", "h/J = 1, run 1 (λ = 1,3,5)", "o"),
    ("zne_time_sweep_ibm_h1_ibm_quebec_lam1-2-3-4-5.json", "h/J = 1, run 2 (λ = 1–5)", "s"),
    ("zne_time_sweep_ibm_h2_ibm_quebec_lam1-2-3-4-5.json", "h/J = 2 (λ = 1–5)", "^"),
]
METHODS = ["raw", "linear", "richardson", "exponential", "exp_fixed"]
LABEL = {"raw": "raw\n(no ZNE)", "linear": "linear ZNE\n(previous\nproject)",
         "richardson": "Richardson", "exponential": "exponential\nfree asym.",
         "exp_fixed": "exponential\nfixed asym."}
SHORT = {"raw": "raw", "linear": "linear ZNE (previous project)", "richardson": "Richardson",
         "exponential": "exponential ZNE", "exp_fixed": "exponential, fixed asymptote"}
# reference palette, categorical slots in fixed order; raw is the neutral baseline
COL = {"raw": "#52514e", "linear": "#eb6834", "richardson": "#1baf7a", "exponential": "#2a78d6",
       "exp_fixed": "#4a3aa7"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#a3a29d", "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titleweight": "bold", "axes.titlesize": 12})


def load():
    runs = []
    for fname, label, marker in RUNS:
        d = json.loads((RES / fname).read_text())
        est = d["observables"]["estimators"]
        per_t = {m: (np.mean([est[m][k]["abs_err"] for k in KEYS], axis=0),
                     np.sqrt(np.sum(np.square([est[m][k]["std"] for k in KEYS]), axis=0)) / len(KEYS))
                 for m in METHODS}
        mean = {m: float(np.mean([est[m][k]["abs_err"] for k in KEYS])) for m in METHODS}
        runs.append({"label": label, "marker": marker, "t": d["t"], "h": d["params"]["h"],
                     "per_t": per_t, "mean": mean})
    return runs


def panel_methods(ax, runs):
    x = np.arange(len(METHODS))
    offsets = np.linspace(-0.22, 0.22, len(runs))
    for r, dx in zip(runs, offsets):
        for i, m in enumerate(METHODS):
            ax.plot(i + dx, r["mean"][m], r["marker"], ms=9, color=COL[m], mec="white", mew=1.5,
                    zorder=3)
    for i, m in enumerate(METHODS):
        vals = [r["mean"][m] for r in runs]
        ax.text(i, max(vals) + 0.006, f"{np.mean(vals):.3f}", ha="center", va="bottom",
                fontsize=10, color=INK, fontweight="bold" if m.startswith("exp") else "normal")
    for r in runs:  # legend: marker shape = run (color carries the method)
        ax.plot([], [], r["marker"], color=MUTED, ms=8, label=r["label"])
    ax.set_xticks(x, [LABEL[m] for m in METHODS], fontsize=9.5)
    ax.set_ylabel("mean |estimate − exact|")
    ax.set_ylim(0, 0.165)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title("(a) Mean error over Mz, Mx, Mzz and t = 0.5–2.0\n     (number = average of the three runs)",
                 loc="left")


def panel_time(ax, run, title, ylabel=True):
    for m in ["raw", "linear", "exponential", "exp_fixed"]:
        mean, std = run["per_t"][m]
        ax.errorbar(run["t"], mean, yerr=std, color=COL[m], marker="o", ms=6, lw=2, capsize=3,
                    label=SHORT[m])
    ax.axvline(1.5, color="#a3a29d", ls=":", lw=1.2)
    ax.text(1.47, 0.205, "previous linear ZNE\nbroke down after t ≈ 1.5", fontsize=8.5, color=MUTED,
            va="top", ha="right")
    ax.set_xlabel("evolution time t (1/J)")
    if ylabel:
        ax.set_ylabel("|estimate − exact|, mean over observables")
    ax.set_ylim(0, 0.21)
    ax.set_xticks(run["t"])
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left")


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    runs = load()
    h1, h2 = runs[1], runs[2]
    t1 = "(b) h/J = 1: error vs time (run 2, λ = 1–5)"
    t2 = "(c) h/J = 2: error vs time (λ = 1–5)"
    foot = ("ibm_quebec, N = 12 periodic ring (0 SWAPs), 2nd-order Trotter dt = 0.1, 4000 shots/circuit, "
            "TREX + Pauli twirling + DD; up to 2400 two-qubit gates at λ = 5")

    fig = plt.figure(figsize=(16, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.28)
    panel_methods(fig.add_subplot(gs[0]), runs)
    ax_b = fig.add_subplot(gs[1])
    panel_time(ax_b, h1, t1)
    ax_b.legend(frameon=False, fontsize=9, loc="center right", bbox_to_anchor=(1.0, 0.42))
    panel_time(fig.add_subplot(gs[2], sharey=ax_b), h2, t2, ylabel=False)
    fig.suptitle("Zero-noise extrapolation on IBM hardware: exponential ZNE cuts the error 3.8–5× "
                 "while linear ZNE removes 17–26%", fontsize=13.5, fontweight="bold", x=0.01, ha="left")
    fig.text(0.01, 0.005, foot, fontsize=9, color=MUTED)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.2)
    p = FIGS / "zne_summary.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.6))
    panel_methods(ax, runs)
    fig.text(0.01, 0.01, foot, fontsize=8.5, color=MUTED)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.92, bottom=0.2)
    fig.savefig(FIGS / "zne_summary_methods.png", dpi=200)
    plt.close(fig)

    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    panel_time(a, h1, t1)
    panel_time(b, h2, t2, ylabel=False)
    a.legend(frameon=False, fontsize=9.5, loc="center right", bbox_to_anchor=(1.0, 0.42))
    fig.text(0.01, 0.01, foot, fontsize=8.5, color=MUTED)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.17, wspace=0.08)
    fig.savefig(FIGS / "zne_summary_vs_time.png", dpi=200)
    plt.close(fig)

    print(p)
    for r in runs:
        print(r["label"], {m: round(v, 4) for m, v in r["mean"].items()})


if __name__ == "__main__":
    main()
