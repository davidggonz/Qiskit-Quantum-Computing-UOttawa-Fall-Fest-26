"""Aer noise models.  (Owner: Toto; initial version ported from Boutaina/Toufic's branch)

  - from_backend(backend): simulator with the backend's calibrated noise (fake or real backend).
  - simple_model(p1q, p2q, p_meas): depolarizing + readout model for quick sweeps
    (defaults are the Iceberg paper's H1-2 values: p1q=4e-4, p2q=3e-3, p_meas=3e-3).
"""
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error


def from_backend(backend):
    """Aer simulator carrying the backend's noise (gate errors, T1/T2, readout)."""
    return AerSimulator.from_backend(backend)


def simple_model(p1q=4e-4, p2q=3e-3, p_meas=3e-3):
    """Depolarizing gate noise + symmetric readout noise.

    `rz` is left noiseless: on IBM hardware it is a virtual (frame-change) gate.
    """
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ["sx", "x", "rx", "h"])
    nm.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ["cz", "ecr", "cx", "rzz"])
    nm.add_all_qubit_readout_error(ReadoutError([[1 - p_meas, p_meas], [p_meas, 1 - p_meas]]))
    return nm


def simulator(backend=None, mode="backend"):
    """AerSimulator for `mode` in {"backend", "simple", "ideal"}."""
    if mode == "backend":
        if backend is None:
            raise ValueError("mode='backend' needs a backend")
        return from_backend(backend)
    if mode == "simple":
        return AerSimulator(noise_model=simple_model())
    if mode == "ideal":
        return AerSimulator()
    raise ValueError(f"unknown noise mode {mode!r}")
