"""Multi-product formulas (MPF).  (Owner: David, with Eli on the coefficient maths)

Idea: run the same evolution time t with several Trotter step counts k_1 < k_2 < ... and
combine the results, <O>_MPF(t) = sum_j x_j <O>_{k_j}(t). For a symmetric 2nd-order formula
the error is a series in (t/k)^2, (t/k)^4, ..., so choosing
    sum_j x_j = 1,   sum_j x_j k_j^(-2m) = 0  for m = 1 .. len(ks)-1
cancels the leading error terms (Richardson extrapolation in the step size).

Cost: the deepest circuit is the one with max(ks) steps, so compare MPF against plain
Trotter with k = max(ks). Shot noise is amplified by ||x||_1 (see `noise_amplification`);
for hardware prefer step counts with small ||x||_1, or the well-conditioned coefficients of
Carrera Vazquez et al., Quantum 7, 1067 (2023) / `qiskit-addon-mpf`.

Caveat: the combination is applied to expectation values, so MPF needs one circuit per k_j.
"""
import numpy as np


def richardson_coefficients(ks, order=2):
    """Static MPF coefficients for step counts `ks`.

    order=2 (symmetric formula): error terms in k^-2, k^-4, ...
    order=1: error terms in k^-1, k^-2, ...
    """
    ks = np.asarray(ks, dtype=float)
    if len(set(ks)) != len(ks):
        raise ValueError("step counts must be distinct")
    p = 2 if order == 2 else 1
    A = np.vstack([ks ** (-p * m) for m in range(len(ks))])
    b = np.zeros(len(ks))
    b[0] = 1.0
    return np.linalg.solve(A, b)


def noise_amplification(coefficients):
    """||x||_1: factor by which independent statistical errors grow in the MPF estimate."""
    return float(np.sum(np.abs(coefficients)))


def mpf_expectation(values_by_k, coefficients):
    """Combine expectation values (scalars or arrays over time), one entry per step count."""
    return sum(c * np.asarray(v, float) for c, v in zip(coefficients, values_by_k))


def mpf_observables(n, J, h, t, ks, order=2, periodic=True):
    """Noiseless MPF estimate of Mz/Mx/Mzz at time t (statevector simulation)."""
    from .tfim import trotter_observables
    x = richardson_coefficients(ks, order)
    runs = [trotter_observables(n, J, h, t, k, order, periodic) for k in ks]
    return {key: float(mpf_expectation([r[key] for r in runs], x)) for key in runs[0]}
