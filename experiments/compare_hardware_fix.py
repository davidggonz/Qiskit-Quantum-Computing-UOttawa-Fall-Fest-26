"""Side-by-side plots: original hardware.py (Boutaina's upload) vs the corrected one.

Loads the original from git (default ref: origin/repo-skeleton) and the corrected file
from the working tree, runs both on the same inputs, and writes four figures:

  1_mz_definition.png   signed <Z> (original) vs RMS Mz (report/team) vs exact, over time
  2_trotter_error.png   |Trotter - ED| on the same observables: original step order vs team 2nd order
  3_zne.png             ZNE: ideal simulator + rzz-only folding (original) vs noisy + folding
  4_layout.png          2-qubit gate count on a heavy-hex fake backend: greedy cluster vs ring

Run:  python experiments/compare_hardware_fix.py [--ref origin/repo-skeleton] [--out results/figs/hardware_fix]
Needs no IBM token (FakeMarrakesh for layout, Aer for noise).
"""
import argparse
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import warnings

import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qfest.ed import quench  # noqa: E402

ORIG, FIXED, EXACT = "#eb6834", "#2a78d6", "#0b0b0b"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.8,
                     "lines.linewidth": 2, "font.size": 10})


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_original(ref):
    src = subprocess.run(["git", "-C", str(ROOT), "show", f"{ref}:hardware.py"],
                         capture_output=True, text=True, check=True).stdout
    tmp = pathlib.Path(tempfile.mkdtemp()) / "hardware_original.py"
    tmp.write_text(src)
    return load_module(tmp, "hardware_original")


def signed_z(psi, n):
    """<Σ Z>/N with Qiskit ordering (qubit i = bit i)."""
    idx = np.arange(2 ** n)
    z = sum(1 - 2 * ((idx >> i) & 1) for i in range(n)) / n
    return float(np.sum(np.abs(psi) ** 2 * z))


def fig_mz_definition(orig, fixed, out, n=8, dt=0.1, steps=50):
    from qiskit.quantum_info import Statevector
    ts = np.arange(steps + 1) * dt
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for ax, h in zip(axes, (0.5, 1.0, 2.0)):
        o_tr = [signed_z(Statevector(orig.build_tfim_circuit(n, 1.0, h, dt, k, add_measure=False)).data, n)
                for k in range(steps + 1)]
        o_ex = [orig.exact_tfim_magnetization(n, 1.0, h, dt, k) for k in range(steps + 1)]
        f_tr = [fixed.observables_from_evs({k: float(np.real(Statevector(
            fixed.build_tfim_circuit(n, 1.0, h, dt, s, add_measure=False)).expectation_value(o)))
            for k, o in fixed.build_tfim_observables(n, 1.0, h).items() if k != "energy"})["Mz"]
            for s in range(steps + 1)]
        ex = quench(n, ts, 1.0, h)["Mz"]
        ax.plot(ts, ex, color=EXACT, lw=2, label="exact Mz (RMS), ring  [reference]")
        ax.plot(ts, f_tr, "o", color=FIXED, ms=4, label="corrected: Mz (RMS), 2nd-order ring")
        ax.plot(ts, o_ex, "--", color=ORIG, lw=2, label="original 'exact': signed ⟨Z⟩, open chain")
        ax.plot(ts, o_tr, "s", color=ORIG, ms=3.5, mfc="none", label="original: signed ⟨Z⟩, 1st-order open")
        ax.set_title(f"h/J = {h}")
        ax.set_xlabel("time t (1/J)")
        ax.axhline(0, color="#a3a29d", lw=0.8)
    axes[0].set_ylabel("magnetization")
    axes[0].legend(loc="lower left", fontsize=8, frameon=False)
    fig.suptitle(f"Mz definition: original signed ⟨Z⟩ vs report's RMS Mz (N = {n}, quench from |0…0⟩, dt = {dt})")
    fig.tight_layout()
    p = out / "1_mz_definition.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_trotter_error(orig, fixed, out, n=8, dt=0.1, steps=50):
    """Same observable, same ring, same ED reference; only the step ordering differs.

    The original's step (ZZ layer, then X layer) equals qfest order=1. Because |0...0>
    is a ZZ eigenstate, (X·ZZ)^k|0> = diagonal phase × (ZZ-outer symmetric formula)|0>,
    so for Z-diagonal observables (Mz, Mzz) it is effectively 2nd order too; Mx is not.
    """
    from qfest.tfim import trotter_curve
    ts = np.arange(steps + 1) * dt
    hs = (1.0, 2.0)
    fig, axes = plt.subplots(len(hs), 3, figsize=(13, 6.4), sharex=True)
    summary = {}
    for row, h in enumerate(hs):
        ex = quench(n, ts, 1.0, h)
        tr1 = trotter_curve(n, 1.0, h, dt, steps, order=1)
        tr2 = trotter_curve(n, 1.0, h, dt, steps, order=2)
        summary[h] = {}
        for col, k in enumerate(("Mz", "Mx", "Mzz")):
            ax = axes[row, col]
            e1, e2 = np.abs(tr1[k] - ex[k]), np.abs(tr2[k] - ex[k])
            ax.semilogy(ts[1:], np.maximum(e1[1:], 1e-12), color=ORIG,
                        label="original step: ZZ then X ('1st order')")
            ax.semilogy(ts[1:], np.maximum(e2[1:], 1e-12), color=FIXED,
                        label="team step: X/2 · ZZ · X/2 (2nd order, fused)")
            ax.set_title(f"{k}, h/J = {h}")
            summary[h][k] = {"orig_step_max": float(e1.max()), "team_2nd_max": float(e2.max())}
        axes[row, 0].set_ylabel("|Trotter − ED|")
    for ax in axes[-1]:
        ax.set_xlabel("time t (1/J)")
    axes[0, 0].legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle(f"Trotter error on the SAME observables (N = {n} ring, dt = {dt}): "
                 "the original ordering is fine for Mz/Mzz, worse for Mx")
    fig.tight_layout()
    p = out / "2_trotter_error.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p, summary


def fig_zne(orig, fixed, out, n=12, h=1.0, dt=0.1, steps=10, shots=4096):
    lams = [1, 3, 5]
    # Original pipeline, exactly as its main() calls it: ideal Aer, folding cx/cz/ecr only.
    core_o = orig.build_tfim_circuit(n, 1.0, h, dt, steps, add_measure=False)
    mag = orig.build_tfim_observables(n, 1.0, h)["magnetization"]
    zo = orig.run_with_zne_aer(core_o, mag, shots=shots)
    exact_o = orig.exact_tfim_magnetization(n, 1.0, h, dt, steps)
    # Corrected pipeline: noisy Aer (Marrakesh medians), fold after transpiling.
    obs = fixed.build_tfim_observables(n, 1.0, h)
    est_obs = {k: obs[k] for k in ("Mz2", "Mx", "Mzz")}
    noise = fixed.noise_model_like(fixed.get_fake_backend("marrakesh"))
    core_f = fixed.build_tfim_circuit(n, 1.0, h, dt, steps, add_measure=False)
    zf = fixed.run_with_zne_aer(core_f, est_obs, shots=shots, noise_model=noise)
    exact_f = fixed.exact_tfim_observables(n, 1.0, h, dt * steps)

    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
    ax = axes[0]
    ax.plot(lams, zo["raw"], "o-", color=ORIG, ms=8, label="raw (ideal sim)")
    ax.plot([0], [zo["mitigated"]], "D", color=ORIG, ms=9, mfc="white", mew=2, label="linear extrapolation")
    ax.axhline(exact_o, color=EXACT, ls="--", lw=1.5, label="its exact (signed ⟨Z⟩)")
    ax.set_title("original: signed ⟨Z⟩")
    ax.set_xlabel("noise factor λ")
    ax.set_xlim(-0.5, 5.5)
    ax.legend(fontsize=8, frameon=False)
    marks = {"linear": "D", "richardson": "^", "exponential": "v"}
    for ax, k in zip(axes[1:], ("Mz", "Mx", "Mzz")):
        ax.plot(lams, zf["raw"][k], "o-", color=FIXED, ms=8, label="raw (noisy sim)")
        for m, mk in marks.items():
            ax.plot([0], [zf["mitigated"][m][k]], mk, color=FIXED, ms=9, mfc="white", mew=2, label=m)
        ax.axhline(exact_f[k], color=EXACT, ls="--", lw=1.5, label="exact (ED)")
        ax.set_title(f"corrected: {k}")
        ax.set_xlabel("noise factor λ")
        ax.set_xlim(-0.5, 5.5)
    axes[1].legend(fontsize=8, frameon=False)
    fig.suptitle(f"ZNE, N = {n} ring, h/J = {h}, t = {dt * steps:g}: original is flat (no noise, nothing folded)")
    fig.tight_layout()
    p = out / "3_zne.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p, {"original": zo, "original_exact": exact_o, "corrected": zf, "corrected_exact": exact_f}


def fig_layout(orig, fixed, out, steps=3):
    from qiskit import transpile
    be = fixed.get_fake_backend("marrakesh")
    rows = []
    for n in (6, 12):
        for periodic in (False, True):
            qc = fixed.build_tfim_circuit(n, 1.0, 1.0, 0.1, steps, periodic=periodic, add_measure=False)
            ideal = 2 * sum(1 for i in qc.data if i.operation.name == "rzz")
            lay_o = orig.select_best_qubits(be, n)
            lay_f = fixed.select_best_qubits(be, n, periodic)
            c = {}
            for label, lay in (("original", lay_o), ("corrected", lay_f)):
                t = transpile(qc, be, optimization_level=3, initial_layout=lay, seed_transpiler=42)
                c[label] = int(t.count_ops().get("cz", 0))
            rows.append({"n": n, "periodic": periodic, "ideal_cz": ideal, **c})
    fig, ax = plt.subplots(figsize=(7.5, 4))
    x = np.arange(len(rows))
    ax.bar(x - 0.2, [r["original"] for r in rows], 0.38, color=ORIG, label="original layout (greedy cluster)")
    ax.bar(x + 0.2, [r["corrected"] for r in rows], 0.38, color=FIXED, label="corrected layout (path / ring)")
    for i, r in enumerate(rows):
        ax.plot([i - 0.42, i + 0.42], [r["ideal_cz"]] * 2, color=EXACT, lw=1.5, ls="--")
        for dx, key in ((-0.2, "original"), (0.2, "corrected")):
            ax.text(i + dx, r[key] + 2, str(r[key]), ha="center", fontsize=8, color="#52514e")
    ax.plot([], [], color=EXACT, ls="--", lw=1.5, label="no-SWAP minimum (2 CZ per RZZ)")
    ax.set_xticks(x, [f"N={r['n']}\n{'ring' if r['periodic'] else 'open'}" for r in rows])
    ax.set_ylabel("CZ count after transpiling")
    ax.set_title(f"Layout on FakeMarrakesh (heavy-hex), {steps} Trotter steps\n"
                 "rings shorter than 12 do not exist on heavy-hex, so N=6 ring always needs SWAPs")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    p = out / "4_layout.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="origin/repo-skeleton", help="git ref holding the original hardware.py")
    ap.add_argument("--out", default=str(ROOT / "results" / "figs" / "hardware_fix"))
    args = ap.parse_args()
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    orig = load_original(args.ref)
    fixed = load_module(ROOT / "hardware.py", "hardware_fixed")

    be = fixed.get_fake_backend("marrakesh")
    summary = {"median_1q_error": {"original": orig.backend_properties_dict(be)["median_1q_err"],
                                   "corrected": fixed.backend_properties_dict(be)["median_1q_err"]}}
    try:
        orig.transpile_compare(orig.build_tfim_circuit(4, add_measure=False), be, levels=(1,))
        summary["original_transpile_compare_with_backend"] = "ok"
    except Exception as ex:
        summary["original_transpile_compare_with_backend"] = f"crash: {type(ex).__name__}: {ex}"

    print(fig_mz_definition(orig, fixed, out))
    p, summary["trotter_max_error"] = fig_trotter_error(orig, fixed, out)
    print(p)
    p, summary["zne"] = fig_zne(orig, fixed, out)
    print(p)
    p, summary["layout"] = fig_layout(orig, fixed, out)
    print(p)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: v for k, v in summary.items() if k != "zne"}, indent=2, default=str))


if __name__ == "__main__":
    main()
