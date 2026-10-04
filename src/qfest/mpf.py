"""Multi-product formulas.  (Owner: David, with Eli on the coefficient maths)

TODO:
  1. Pick step counts, e.g. ks = [1, 2, 4] for a total time t (Trotter with dt = t/k).
  2. Solve for MPF coefficients x_j (Richardson-style: sum x_j = 1, sum x_j k_j^{-2m}=0 ...).
     Start with the plain Richardson coefficients; compare with the well-conditioned
     variant (Carrera Vazquez et al., Quantum 7, 1067) and `qiskit-addon-mpf`.
  3. Estimate <O>(t) = sum_j x_j <O>_j(t) from separate circuits (one per k_j).
Acceptance: beats plain 2nd-order Trotter at equal max depth for h/J = 2 (report: 8.46% at dt=0.05).
"""


def richardson_coefficients(ks, order=2):
    raise NotImplementedError


def mpf_expectation(values_by_k, coefficients):
    raise NotImplementedError
