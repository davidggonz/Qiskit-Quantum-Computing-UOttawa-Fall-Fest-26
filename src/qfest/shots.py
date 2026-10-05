"""Finite-shot estimators from statevectors (Z-basis sampling).  Shared helper.

Mzz and Mz^2 are diagonal in the computational basis, so one Z-basis measurement gives both.
Sampling counts from |psi|^2 mimics what a noiseless device with `shots` shots would return.
"""
import numpy as np

from .ed import observables_ops


def z_diagonals(n, periodic=True):
    """Per-basis-state values of the diagonal observables {"Mzz", "Mz2"}."""
    ops = observables_ops(n, periodic)
    return {k: np.real(ops[k].diagonal()) for k in ("Mzz", "Mz2")}


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
