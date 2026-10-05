"""Exact-diagonalization baseline for the periodic/open TFIM.  (Owner: Ririsha)

H = -J sum Z_i Z_{i+1} - h sum X_i.  Qubit i is the i-th least-significant bit,
matching Qiskit's little-endian convention, so statevectors compare directly.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, expm_multiply

_X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex)
_Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)


def site_op(op, i, n):
    """Single-site operator `op` acting on qubit i of n."""
    return sp.kron(sp.kron(sp.identity(2 ** (n - 1 - i)), op), sp.identity(2 ** i), format="csr")


def bonds(n, periodic=True):
    b = [(i, i + 1) for i in range(n - 1)]
    if periodic and n > 2:
        b.append((n - 1, 0))
    return b


def hamiltonian(n, J=1.0, h=1.0, periodic=True):
    H = sp.csr_matrix((2 ** n, 2 ** n), dtype=complex)
    for i, j in bonds(n, periodic):
        H = H - J * site_op(_Z, i, n) @ site_op(_Z, j, n)
    for i in range(n):
        H = H - h * site_op(_X, i, n)
    return H


def observables_ops(n, periodic=True):
    """Return dict of sparse operators: Mz2 (=(sum Z)^2 / N^2), Mx, Mzz."""
    Zs = sum((site_op(_Z, i, n) for i in range(n)), sp.csr_matrix((2 ** n,) * 2, dtype=complex))
    Mzz = sum((site_op(_Z, i, n) @ site_op(_Z, j, n) for i, j in bonds(n, periodic)),
              sp.csr_matrix((2 ** n,) * 2, dtype=complex))
    Mx = sum((site_op(_X, i, n) for i in range(n)), sp.csr_matrix((2 ** n,) * 2, dtype=complex))
    return {"Mz2": (Zs @ Zs) / n ** 2, "Mx": Mx / n, "Mzz": Mzz / n}


def measure(psi, ops):
    """Observables as in the report: Mz is the RMS magnitude sqrt(<(sum Z)^2>)/N."""
    ev = {k: float(np.real(np.vdot(psi, O @ psi))) for k, O in ops.items()}
    return {"Mz": float(np.sqrt(max(ev["Mz2"], 0.0))), "Mx": ev["Mx"], "Mzz": ev["Mzz"]}


def ground_state(n, J=1.0, h=1.0, periodic=True):
    H = hamiltonian(n, J, h, periodic)
    w, v = eigsh(H, k=1, which="SA")
    return w[0].real, v[:, 0]


def ground_state_observables(n, J=1.0, h=1.0, periodic=True):
    _, psi = ground_state(n, J, h, periodic)
    return measure(psi, observables_ops(n, periodic))


def quench(n, times, J=1.0, h=1.0, periodic=True, psi0=None):
    """Exact evolution from |0...0> (default). Returns dict of arrays over `times`."""
    H = hamiltonian(n, J, h, periodic)
    ops = observables_ops(n, periodic)
    if psi0 is None:
        psi0 = np.zeros(2 ** n, dtype=complex)
        psi0[0] = 1.0
    if len(times) == 1:
        states = [expm_multiply(-1j * float(times[0]) * H, psi0)]
    else:  # times must be uniformly spaced
        states = expm_multiply(-1j * H, psi0, start=float(times[0]), stop=float(times[-1]),
                               num=len(times), endpoint=True)
    out = {"t": np.asarray(times, dtype=float), "Mz": [], "Mx": [], "Mzz": []}
    for psi in states:
        m = measure(psi, ops)
        for k in ("Mz", "Mx", "Mzz"):
            out[k].append(m[k])
    return {k: np.asarray(v) for k, v in out.items()}
