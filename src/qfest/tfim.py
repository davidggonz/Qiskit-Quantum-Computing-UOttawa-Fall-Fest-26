"""TFIM Trotter circuit builder.  (Owner: David)

Interface contract (do not change signature without telling the team):
    tfim_circuit(n, J, h, dt, steps, order=2, periodic=True) -> QuantumCircuit
Initial state is |0...0> (quench protocol).  No measurements are added.
Order 2 fuses adjacent half-step Rx rotations (r+1 Rx layers instead of 2r).
"""
from qiskit import QuantumCircuit
from .ed import bonds


def edge_colouring(n, periodic=True):
    """Group bonds into layers of disjoint edges (2 layers for even rings/chains, 3 for odd rings)."""
    layers = {}
    for k, (i, j) in enumerate(bonds(n, periodic)):
        c = k % 2
        if periodic and n % 2 == 1 and k == n - 1:
            c = 2
        layers.setdefault(c, []).append((i, j))
    return [layers[c] for c in sorted(layers)]


def _zz_layer(qc, layers, J, dt):
    for layer in layers:
        for i, j in layer:
            qc.rzz(-2.0 * J * dt, i, j)  # exp(+i J dt Z Z)


def _rx_layer(qc, n, h, angle_dt):
    for i in range(n):
        qc.rx(-2.0 * h * angle_dt, i)  # exp(+i h angle_dt X)


def tfim_circuit(n, J=1.0, h=1.0, dt=0.05, steps=1, order=2, periodic=True):
    qc = QuantumCircuit(n)
    layers = edge_colouring(n, periodic)
    if order == 1:
        for _ in range(steps):
            _zz_layer(qc, layers, J, dt)
            _rx_layer(qc, n, h, dt)
    elif order == 2:
        _rx_layer(qc, n, h, dt / 2)
        for s in range(steps):
            _zz_layer(qc, layers, J, dt)
            _rx_layer(qc, n, h, dt if s < steps - 1 else dt / 2)
    else:
        raise ValueError("order must be 1 or 2")
    return qc


def _step_unitary(n, J, h, dt, order, periodic):
    from qiskit.quantum_info import Operator
    return Operator(tfim_circuit(n, J, h, dt, 1, order, periodic)).data


def trotter_curve(n, J=1.0, h=1.0, dt=0.05, steps=1, order=2, periodic=True):
    """Noiseless Trotter observables at t = 0, dt, ..., steps*dt (statevector, small n only).

    Repeats the single-step unitary, which equals the fused circuit at step boundaries.
    Returns {"t", "Mz", "Mx", "Mzz"} as numpy arrays.
    """
    import numpy as np
    from .ed import observables_ops, measure
    ops = observables_ops(n, periodic)
    U = _step_unitary(n, J, h, dt, order, periodic)
    psi = np.zeros(2 ** n, dtype=complex)
    psi[0] = 1.0
    out = {"Mz": [], "Mx": [], "Mzz": []}
    for s in range(steps + 1):
        m = measure(psi, ops)
        for k in out:
            out[k].append(m[k])
        if s < steps:
            psi = U @ psi
    out = {k: np.asarray(v) for k, v in out.items()}
    out["t"] = np.arange(steps + 1) * dt
    return out


def trotter_observables(n, J=1.0, h=1.0, t=1.0, steps=1, order=2, periodic=True):
    """Noiseless observables after `steps` Trotter steps of size t/steps (statevector)."""
    from qiskit.quantum_info import Statevector
    from .ed import observables_ops, measure
    psi = Statevector(tfim_circuit(n, J, h, t / steps, steps, order, periodic)).data
    return measure(psi, observables_ops(n, periodic))
