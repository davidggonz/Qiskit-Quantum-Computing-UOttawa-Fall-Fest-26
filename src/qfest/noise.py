"""Aer noise models.  (Owner: Toto)

TODO:
  - from_backend(backend): AerSimulator/NoiseModel.from_backend from a fake or real backend.
  - simple_model(p1q, p2q, p_meas): depolarizing + readout model for quick sweeps
    (the Iceberg paper's H1-2 values: p1q=4e-4, p2q=3e-3, p_meas=3e-3).
"""


def from_backend(backend):
    raise NotImplementedError


def simple_model(p1q=4e-4, p2q=3e-3, p_meas=3e-3):
    raise NotImplementedError
