"""Figures for the README results section, built from the saved IBM hardware runs.

Writes to results/figs/:
  readme_error_by_method.png   mean error per run: raw vs linear ZNE vs exponential ZNE
  readme_observables_vs_time.png   Mz, Mx, Mzz vs time at h/J = 1 against the exact curve
  readme_validation.png        run-to-run reproducibility and chi-squared fit tests

The exact curve is an independent dense exact diagonalization of the N = 12 ring; the
script checks it against the exact values stored with the hardware data.

Run:  python experiments/plot_readme_figures.py
"""
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES, FIGS = ROOT / "results", ROOT / "results" / "figs"
KEYS = ["Mz", "Mx", "Mzz"]
RUNS = [  # (file, label)
    ("zne_time_sweep_ibm_h1_ibm_quebec.json", "h/J = 1, run 1"),
    ("zne_time_sweep_ibm_h1_ibm_quebec_lam1-2-3-4-5.json", "h/J = 1, run 2"),
    ("zne_time_sweep_ibm_h2_ibm_quebec_lam1-2-3-4-5.json", "h/J = 2"),
]
COL = {"exact": "#0b0b0b", "raw": "#8a8984", "linear": "#eb6834", "exponential": "#2a78d6"}
NAME = {"raw": "Raw (no mitigation)", "linear": "Linear ZNE (previous method)",
        "exponential": "Exponential ZNE (this work)"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CHI2_CUT = 5.991 / 2  # chi^2/dof at the 5% significance level for 2 degrees of freedom

plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#a3a29d", "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titleweight": "bold", "axes.titlesize": 13})


def load(fname):
    return json.loads((RES / fname).read_text())


def mean_err(d, m):
    return float(np.mean([d["observables"]["estimators"][m][k]["abs_err"] for k in KEYS]))


def exact_curve(n, J, h, times):
    """Dense ED of H = -J sum Z_i Z_{i+1} - h sum X_i (periodic), quench from |0...0>."""
    dim = 2 ** n
    idx = np.arange(dim)
    z = 1 - 2 * ((idx[:, None] >> np.arange(n)[None, :]) & 1)  # z[s, i] = +-1
    zz = (z * np.roll(z, -1, axis=1)).sum(axis=1)
    H = np.diag(-J * zz.astype(float))
    for i in range(n):
        H[idx, idx ^ (1 << i)] += -h
    w, v = np.linalg.eigh(H)
    c0 = v[0, :]  # overlaps <E_n|0...0>
    sz = z.sum(axis=1).astype(float)
    out = {k: [] for k in KEYS}
    for t in times:
        psi = v @ (np.exp(-1j * w * t) * c0)
        p = np.abs(psi) ** 2
        out["Mz"].append(np.sqrt(p @ sz ** 2) / n)
        out["Mzz"].append(p @ zz / n)
        xpsi = np.zeros_like(psi)
        for i in range(n):
            xpsi += psi[idx ^ (1 << i)]
        out["Mx"].append(np.real(np.vdot(psi, xpsi)) / n)
    return {k: np.array(v_) for k, v_ in out.items()}


def fig_error_by_method(runs):
    fig, ax = plt.subplots(figsize=(11, 5.6))
    methods = ["raw", "linear", "exponential"]
    x = np.arange(len(runs))
    width = 0.26
    for j, m in enumerate(methods):
        vals = [mean_err(d, m) for d, _ in runs]
        bars = ax.bar(x + (j - 1) * width, vals, width, color=COL[m], label=NAME[m])
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.003, f"{v:.3f}", ha="center",
                    va="bottom", fontsize=11, color=INK,
                    fontweight="bold" if m == "exponential" else "normal")
    for i, (d, _) in enumerate(runs):
        r = mean_err(d, "raw") / mean_err(d, "exponential")
        ax.annotate(f"{r:.1f}x smaller", xy=(x[i] + width, mean_err(d, "exponential") + 0.03),
                    ha="center", fontsize=11, color=COL["exponential"], fontweight="bold")
    ax.set_xticks(x, [lab for _, lab in runs])
    ax.set_ylabel("mean |estimate - exact|  (lower is better)")
    ax.set_ylim(0, 0.175)
    ax.grid(axis="y", color=GRID)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("Error on ibm_quebec, averaged over Mz, Mx, Mzz and t = 0.5-2.0", loc="left")
    fig.tight_layout()
    fig.savefig(FIGS / "readme_error_by_method.png", dpi=200)
    plt.close(fig)


def fig_observables(d):
    p = d["params"]
    tt = np.linspace(0, 2.1, 85)
    ex = exact_curve(p["n"], p["J"], p["h"], tt)
    chk = exact_curve(p["n"], p["J"], p["h"], d["t"])
    for k in KEYS:  # independent ED must reproduce the stored exact values
        assert np.allclose(chk[k], d["observables"]["exact"][k], atol=1e-6), k
    est = d["observables"]["estimators"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.9), sharey=True)
    labels = {"Mz": "Mz (RMS magnetization)", "Mx": "Mx (transverse magnetization)",
              "Mzz": "Mzz (nearest-neighbour correlation)"}
    for ax, k in zip(axes, KEYS):
        ax.plot(tt, ex[k], color=COL["exact"], lw=2, label="Exact (diagonalization)")
        for m, off, mk in [("raw", -0.03, "o"), ("linear", 0.0, "s"), ("exponential", 0.03, "D")]:
            ax.errorbar(np.array(d["t"]) + off, est[m][k]["mean"], yerr=est[m][k]["std"],
                        fmt=mk, ms=7, color=COL[m], capsize=3, lw=1.5, label=NAME[m])
        ax.set_title(labels[k], loc="left")
        ax.set_xlabel("time t (1/J)")
        ax.set_xticks([0, 0.5, 1.0, 1.5, 2.0])
        ax.grid(color=GRID)
        ax.set_axisbelow(True)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("value")
    axes[0].legend(frameon=False, fontsize=10, loc="lower left")
    fig.suptitle("h/J = 1 on ibm_quebec: exponential ZNE follows the exact dynamics, "
                 "linear ZNE stays near the raw data", x=0.01, ha="left", fontsize=13.5,
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGS / "readme_observables_vs_time.png", dpi=200)
    plt.close(fig)


def raw_lambda1(d):
    """Raw (lambda = 1) Mz, Mx, Mzz per time; raw_evs rows are (Mz^2, Mx, Mzz)."""
    out = []
    for row, ev in zip(d["extra"]["index"], d["extra"]["raw_evs"]):
        if float(row["lambda"]) == 1.0:
            out.extend([np.sqrt(ev[0]), ev[1], ev[2]])
    return np.array(out)


def fig_validation(runs):
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1, 1.5]})
    r1, r2 = raw_lambda1(runs[0][0]), raw_lambda1(runs[1][0])
    a.plot([0, 0.8], [0, 0.8], color=MUTED, lw=1, ls="--", label="perfect agreement")
    for j, k in enumerate(KEYS):
        a.plot(r1[j::3], r2[j::3], "o", ms=8, label=k)
    a.set_xlabel("run 1 (raw value)")
    a.set_ylabel("run 2, a day later (raw value)")
    a.set_title(f"Reproducible: max difference {np.max(np.abs(r1 - r2)):.3f}", loc="left")
    a.set_xlim(0, 0.8)
    a.set_ylim(0, 0.8)
    a.grid(color=GRID)
    a.legend(frameon=False, fontsize=10)

    names, vals = [], []
    for d, lab in runs[1:]:
        ch = d["extra"]["chi2_and_dof"]["exponential"]
        for i, t in enumerate(d["t"]):
            for k in KEYS:
                c, dof = ch[k][i]
                names.append(f"{k} t={t}")
                vals.append(c / dof)
    vals = np.array(vals)
    ok = vals <= CHI2_CUT
    xs = np.arange(len(vals))
    b.bar(xs, np.minimum(vals, 9), color=np.where(ok, COL["exponential"], "#c9c8c3"))
    b.axhline(CHI2_CUT, color=COL["linear"], ls="--", lw=1.5,
              label=f"5% significance limit (chi2/dof = {CHI2_CUT:.2f})")
    b.axhline(1, color=MUTED, ls=":", lw=1, label="ideal fit (chi2/dof = 1)")
    b.axvline(11.5, color="#a3a29d", lw=1)
    b.text(5.5, 8.4, "h/J = 1", ha="center", color=MUTED)
    b.text(17.5, 8.4, "h/J = 2", ha="center", color=MUTED)
    b.set_xticks(xs, names, rotation=90, fontsize=8)
    b.set_ylabel("chi2 per degree of freedom (capped at 9)")
    b.set_ylim(0, 9)
    b.set_title(f"Exponential model passes the fit test at {ok.sum()} of {len(vals)} points",
                loc="left")
    b.legend(frameon=False, fontsize=10, loc="upper right", bbox_to_anchor=(1.0, 0.93))
    fig.tight_layout()
    fig.savefig(FIGS / "readme_validation.png", dpi=200)
    plt.close(fig)
    return ok.sum(), len(vals), float(np.max(np.abs(r1 - r2)))


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    runs = [(load(f), lab) for f, lab in RUNS]
    fig_error_by_method(runs)
    fig_observables(runs[1][0])
    print("validation:", fig_validation(runs))
    for d, lab in runs:
        print(lab, {m: round(mean_err(d, m), 3) for m in ["raw", "linear", "exponential"]})


if __name__ == "__main__":
    main()
