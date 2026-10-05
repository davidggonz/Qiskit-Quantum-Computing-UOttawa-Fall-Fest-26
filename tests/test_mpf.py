import numpy as np
import pytest

from qfest.ed import quench
from qfest.mpf import richardson_coefficients, mpf_observables, noise_amplification
from qfest.tfim import trotter_observables


def test_coefficients_cancel_error_terms():
    ks = [2, 3, 4]
    x = richardson_coefficients(ks, order=2)
    assert x.sum() == pytest.approx(1.0)
    assert np.dot(x, np.array(ks, float) ** -2) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(x, np.array(ks, float) ** -4) == pytest.approx(0.0, abs=1e-12)
    assert noise_amplification(x) >= 1.0


def test_mpf_beats_trotter_at_same_max_depth():
    n, h, t, ks = 6, 2.0, 1.0, [4, 6, 8]
    exact = quench(n, [t], 1.0, h)["Mzz"][0]
    trotter_err = abs(trotter_observables(n, 1.0, h, t, max(ks))["Mzz"] - exact)
    mpf_err = abs(mpf_observables(n, 1.0, h, t, ks)["Mzz"] - exact)
    assert mpf_err < trotter_err / 10


def test_table2_reproduced_exactly():
    """Old report Table 2 entries for h/J = 1 (the most sensitive row)."""
    from qfest.metrics import max_report_deviation
    from qfest.tfim import trotter_curve
    for dt, expected in ((0.2, 50.53), (0.1, 12.66), (0.05, 3.12)):
        tr = trotter_curve(6, 1.0, 1.0, dt, int(round(20 / dt)))
        ex = quench(6, tr["t"], 1.0, 1.0)
        assert max_report_deviation(tr["Mzz"][1:], ex["Mzz"][1:]) == pytest.approx(expected, abs=0.006)
