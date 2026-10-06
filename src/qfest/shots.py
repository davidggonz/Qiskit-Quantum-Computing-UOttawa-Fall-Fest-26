"""Finite-shot estimators from statevectors (Z-basis sampling).  Shared helper.

Mzz and Mz^2 are diagonal in the computational basis, so one Z-basis measurement gives both.
Sampling counts from |psi|^2 mimics what a noiseless device with `shots` shots would return.
"""
import numpy as np

from .ed import observables_ops


def z_diagonals(n, periodic=True):
    """Per-basis-state values of the diagonal observables {"Mzz", "Mz2", "Zavg"}.

    Zavg = (1/N) sum_i Z_i is the signed magnetization; Mz2 = (sum_i Z_i)^2 / N^2, so the
    report's RMS Mz is sqrt(<Mz2>). Basis index bit i = qubit i (Qiskit order).
    """
    ops = observables_ops(n, periodic)
    out = {k: np.real(ops[k].diagonal()) for k in ("Mzz", "Mz2")}
    idx = np.arange(2 ** n)
    out["Zavg"] = sum(1.0 - 2.0 * ((idx >> i) & 1) for i in range(n)) / n
    return out


def counts_expectations(counts, n, periodic=True, diags=None):
    """Zavg, Mzz and RMS Mz from a Qiskit counts dict whose classical bit i holds qubit i."""
    diags = diags or z_diagonals(n, periodic)
    shots = sum(counts.values())
    ev = {k: 0.0 for k in diags}
    for bitstring, c in counts.items():
        j = int(bitstring.replace(" ", ""), 2)
        for k, d in diags.items():
            ev[k] += c * d[j]
    ev = {k: v / shots for k, v in ev.items()}
    return {"Zavg": ev["Zavg"], "Mzz": ev["Mzz"], "Mz": float(np.sqrt(max(ev["Mz2"], 0.0)))}


def state_expectations(psi, n, periodic=True, diags=None):
    """Same observables as `counts_expectations`, exactly, from a statevector."""
    diags = diags or z_diagonals(n, periodic)
    p = np.abs(np.asarray(psi)) ** 2
    ev = {k: float(p @ d) for k, d in diags.items()}
    return {"Zavg": ev["Zavg"], "Mzz": ev["Mzz"], "Mz": float(np.sqrt(max(ev["Mz2"], 0.0)))}


def sample_estimate(psi, diag, shots, rng):
    """Shot estimate of a diagonal observable given its per-basis values `diag`."""
    p = np.abs(psi) ** 2
    counts = rng.multinomial(shots, p / p.sum())
    return float(counts @ diag) / shots


def allocate_shots(coefficients, total_shots):
    """Split a shot budget across MPF circuits proportional to |x_j| (minimises variance)."""
    w = np.abs(np.asarray(coefficients, float))
    s = np.maximum(1, np.round(total_shots * w / w.sum()).astype(int))
    return s
