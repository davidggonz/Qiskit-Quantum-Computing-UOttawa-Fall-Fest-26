import numpy as np
import pytest

from qfest.shots import z_diagonals, sample_estimate, allocate_shots
from qfest.tfim import trotter_state


def test_shot_estimate_is_unbiased():
    psi = trotter_state(4, 1.0, 1.0, 0.7, 4)
    diag = z_diagonals(4)["Mzz"]
    exact = float(np.real(np.vdot(psi, diag * psi)))
    rng = np.random.default_rng(0)
    est = np.mean([sample_estimate(psi, diag, 2000, rng) for _ in range(200)])
    assert est == pytest.approx(exact, abs=5e-3)


def test_allocation_follows_weights():
    s = allocate_shots([2.0, -1.0, 0.0], 3000)
    assert s[0] == pytest.approx(2 * s[1], rel=0.01) and s[2] >= 1
