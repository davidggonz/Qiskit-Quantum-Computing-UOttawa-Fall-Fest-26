"""Re-run Boutaina's canonical setup with the ported hardware module and plot both side by side.

Her canonical run (feature/canonical-tfim, results_canonical/tfim_results.json):
  N=6 periodic ring, h/J=0.5, dt=0.05, 20 steps, quench from |0..0>, FakeMarrakesh noise,
  observable = signed <Z> (her "M_z"), ZNE scales 1,3,5 with linear extrapolation.
The port cannot place a 6-site ring without SWAPs on heavy-hex, so by default it runs the
6-site open chain (same h, dt, steps) and, with --ring12, also the 12-site ring.

Produces results/figs/compare_{accuracy,depth,zne}.png and results/compare_hardware_runs.json.

  python experiments/compare_hardware_runs.py --boutaina path/to/results_canonical/tfim_results.json
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qfest import hardware as hw
from qfest.ed import quench_state
from qfest.extrapolation import EXTRAPOLATORS
from qfest.noise import simulator
from qfest.results import RESULTS_DIR, save
from qfest.shots import counts_expectations, state_expectations
from qfest.tfim import tfim_circuit, trotter_state

# validated categorical palette (dataviz reference instance), fixed order
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
LEVELS = (0, 1, 2, 3)


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    ax.set_xlabel(xlabel, color=INK2)
    ax.set_ylabel(ylabel, color=INK2)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2)


def run_port(n, periodic, h, dt, steps, shots, seed, fake):
    backend = hw.get_backend(fake=fake)
    layout = hw.find_layout(backend, n, periodic)
    logical = tfim_circuit(n, 1.0, h, dt, steps, periodic=periodic)
    sim = simulator(backend, "backend")
    t = dt * steps
    exact = state_expectations(quench_state(n, t, 1.0, h, periodic), n, periodic)
    ideal = state_expectations(trotter_state(n, 1.0, h, t, steps, periodic=periodic), n, periodic)
    levels = []
    for lvl in LEVELS:
        isa = hw.to_isa(logical, backend, layout, lvl, seed)
        rep = hw.transpile_report(logical, isa)
        circ = hw.add_measurements(isa, hw.physical_qubits(isa)[:n])
        val = counts_expectations(hw.run_simulator([circ], sim, shots, seed + lvl)[0], n, periodic)["Zavg"]
        levels.append({"optimization_level": lvl, "depth": rep["depth"], "two_qubit": rep["two_qubit"],
                       "extra_two_qubit": rep["extra_two_qubit"], "Zavg": val, "abs_error": abs(val - exact["Zavg"])})
    isa = hw.to_isa(logical, backend, layout, 1, seed)
    circs, factors = hw.zne_circuits(isa, n, backend, (1, 3, 5))
    raw = [counts_expectations(c, n, periodic)["Zavg"] for c in hw.run_simulator(circs, sim, shots, seed)]
    return {"label": f"Port: N={n} {'ring' if periodic else 'chain'}, zero-SWAP layout",
            "n": n, "periodic": periodic, "layout": layout, "exact": exact["Zavg"], "ideal": ideal["Zavg"],
            "levels": levels, "zne": {"factors": factors, "raw": raw,
                                      "mitigated": {k: f(factors, raw) for k, f in EXTRAPOLATORS.items()}}}


def load_boutaina(path):
    d = json.loads(Path(path).read_text())
    rows = [r for r in d["transpilation"]]
    aer = d["aer_expectations"]
    exact = d["exact_magnetization"]
    raw = d["zne_aer"]["raw"]
    factors = [float(x) for x in d["zne_aer"]["noise_factors"]]
    return {"label": "Boutaina: N=6 ring, auto layout (SWAPs)", "n": 6, "periodic": True,
            "exact": exact, "ideal": None,
            "levels": [{"optimization_level": r["optimization_level"], "depth": r["depth"], "two_qubit": r["cx_like"],
                        "Zavg": aer[str(r["optimization_level"])],
                        "abs_error": abs(aer[str(r["optimization_level"])] - exact)} for r in rows],
            "zne": {"factors": factors, "raw": raw,
                    "mitigated": {k: f(factors, raw) for k, f in EXTRAPOLATORS.items()}}}


def plot_accuracy(runs, out):
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    for run, color, marker in zip(runs, (ORANGE, BLUE, AQUA), ("o", "s", "D")):
        x = [r["optimization_level"] for r in run["levels"]]
        y = [r["abs_error"] for r in run["levels"]]
        ax.plot(x, y, color=color, linewidth=2, marker=marker, markersize=8, label=run["label"])
        ax.annotate(f"{y[-1]:.3f}", (x[-1], y[-1]), textcoords="offset points", xytext=(8, 0),
                    va="center", fontsize=8, color=INK2)
    ax.set_xticks(LEVELS)
    ax.set_ylim(bottom=0)
    style(ax, "Accuracy vs optimization level (signed ⟨Z⟩, FakeMarrakesh noise)",
          "optimization level", "|noisy − exact|")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(out, facecolor="#fcfcfb")
    plt.close(fig)


def plot_depth(runs, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), dpi=150)
    width = 0.8 / len(runs)
    for k, (key, title) in enumerate((("depth", "Transpiled depth"), ("two_qubit", "Two-qubit gates (CZ)"))):
        ax = axes[k]
        for i, (run, color) in enumerate(zip(runs, (ORANGE, BLUE, AQUA))):
            x = np.array(LEVELS) + (i - (len(runs) - 1) / 2) * width
            y = [r[key] for r in run["levels"]]
            bars = ax.bar(x, y, width * 0.92, color=color, label=run["label"], edgecolor="#fcfcfb", linewidth=1)
            for b, v in zip(bars, y):
                ax.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=7, color=INK2)
        ax.set_xticks(LEVELS)
        style(ax, title, "optimization level", key.replace("_", "-"))
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left", bbox_to_anchor=(0, -0.18), ncol=1)
    fig.tight_layout()
    fig.savefig(out, facecolor="#fcfcfb", bbox_inches="tight")
    plt.close(fig)


def plot_zne(runs, out):
    fig, axes = plt.subplots(1, len(runs), figsize=(4.8 * len(runs), 4.4), dpi=150, sharey=True)
    axes = np.atleast_1d(axes)
    ex_colors = {"linear": ORANGE, "richardson": AQUA, "exponential": YELLOW}
    for ax, run in zip(axes, runs):
        z = run["zne"]
        ax.axhline(run["exact"], color=INK2, linewidth=1.5, linestyle="--")
        ax.text(5.4, run["exact"], f"exact {run['exact']:.3f}", va="bottom", ha="right", fontsize=8, color=INK2)
        ax.plot(z["factors"], z["raw"], color=BLUE, linewidth=2, marker="o", markersize=8, label="raw (noisy)")
        for name, val in z["mitigated"].items():
            ax.plot([0], [val], marker="*", markersize=14, color=ex_colors[name], linestyle="none",
                    markeredgecolor="#fcfcfb", label=f"{name}: {val:.3f}")
        ax.set_xlim(-0.4, 5.6)
        ax.set_xticks([0, 1, 3, 5])
        style(ax, run["label"], "noise factor", "signed ⟨Z⟩")
        ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower left")
    axes[0].set_ylim(0, 1.0)
    fig.tight_layout()
    fig.savefig(out, facecolor="#fcfcfb")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--boutaina", default=None, help="her results_canonical/tfim_results.json")
    p.add_argument("--ring12", action="store_true", help="also run the 12-site ring (~5 min)")
    p.add_argument("--h", type=float, default=0.5)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--steps", type=int, default=20)
    p.add_argument("--shots", type=int, default=4096)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--fake", default="marrakesh")
    a = p.parse_args()

    runs = [load_boutaina(a.boutaina)] if a.boutaina else []
    runs.append(run_port(6, False, a.h, a.dt, a.steps, a.shots, a.seed, a.fake))
    if a.ring12:
        runs.append(run_port(12, True, a.h, a.dt, a.steps, a.shots, a.seed, a.fake))

    figs = RESULTS_DIR / "figs"
    figs.mkdir(parents=True, exist_ok=True)
    plot_accuracy(runs, figs / "compare_accuracy.png")
    plot_depth(runs, figs / "compare_depth.png")
    plot_zne(runs, figs / "compare_zne.png")
    save({"experiment": "compare_hardware_runs", "backend": f"aer_backend_fake_{a.fake}",
          "params": vars(a), "method": "trotter2_zne", "observables": {}, "extra": {"runs": runs}},
         "compare_hardware_runs")
    for run in runs:
        z = run["zne"]
        print(f"{run['label']}: exact {run['exact']:.4f}, raw x1 {z['raw'][0]:.4f}, "
              + ", ".join(f"{k} {v:.4f}" for k, v in z["mitigated"].items()))
    print(f"figures in {figs}")


if __name__ == "__main__":
    main()
