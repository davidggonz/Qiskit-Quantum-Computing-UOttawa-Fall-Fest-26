import numpy as np
import pytest

from qfest.ed import ground_state_observables, quench
from qfest.extrapolation import linear, richardson, exponential

# Entries of the previous report's Table 1 that are self-consistent (periodic chain, J=1).
# Other entries of that table disagree with exact diagonalization and with Kramers-Wannier
# duality (Mx(h) = Mzz(1/h)); see docs/TABLE1_DISCREPANCIES.md.
TABLE1_OK = [
    (6, 0.5, "Mz", 0.969914), (6, 0.5, "Mzz", 0.931517),
    (6, 2.0, "Mx", 0.931517), (6, 2.0, "Mzz", 0.265198),
    (8, 0.5, "Mx", 0.260062), (8, 0.5, "Mzz", 0.933604),
    (8, 2.0, "Mz", 0.481939), (8, 2.0, "Mx", 0.933604), (8, 2.0, "Mzz", 0.260062),
    (10, 2.0, "Mx", 0.934076), (10, 2.0, "Mzz", 0.258970),
]


@pytest.mark.parametrize("n,h,key,val", TABLE1_OK)
def test_ground_state_matches_report(n, h, key, val):
    assert ground_state_observables(n, 1.0, h)[key] == pytest.approx(val, abs=2e-6)


@pytest.mark.parametrize("n", [6, 8])
def test_kramers_wannier_duality(n):
    lo, hi = ground_state_observables(n, 1.0, 0.5), ground_state_observables(n, 1.0, 2.0)
    assert lo["Mx"] == pytest.approx(hi["Mzz"], abs=1e-8)
    assert lo["Mzz"] == pytest.approx(hi["Mx"], abs=1e-8)
    crit = ground_state_observables(n, 1.0, 1.0)
    assert crit["Mx"] == pytest.approx(crit["Mzz"], abs=1e-8)


def test_trotter_circuit_converges_to_ed():
    from qiskit.quantum_info import Statevector
    from qfest.tfim import tfim_circuit
    from qfest.ed import observables_ops, measure
    n, h, T = 6, 0.5, 2.0
    ops = observables_ops(n)
    exact = quench(n, [T], 1.0, h)["Mzz"][0]
    errs = []
    for dt in (0.2, 0.1, 0.05):
        psi = Statevector(tfim_circuit(n, 1.0, h, dt, int(round(T / dt)), order=2)).data
        errs.append(abs(measure(psi, ops)["Mzz"] - exact))
    assert errs[0] > errs[1] > errs[2]
    assert errs[2] < 1e-3


def test_extrapolators_recover_known_decay():
    lams = [1, 3, 5]
    vals = [0.8 * np.exp(-0.1 * x) for x in lams]
    assert exponential(lams, vals, asymptote=0.0) == pytest.approx(0.8, abs=1e-6)
    assert abs(linear(lams, vals) - 0.8) > abs(richardson(lams, vals) - 0.8)
