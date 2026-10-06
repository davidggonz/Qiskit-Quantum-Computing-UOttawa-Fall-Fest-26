"""Checks for the ported hardware module on an offline fake Heron backend."""
import numpy as np
import pytest

pytest.importorskip("qiskit_ibm_runtime")

from qfest import hardware as hw
from qfest.noise import simulator
from qfest.shots import counts_expectations, state_expectations
from qfest.tfim import tfim_circuit, trotter_state


@pytest.fixture(scope="module")
def backend():
    return hw.get_backend(fake="marrakesh")


@pytest.mark.parametrize("n,periodic", [(12, True), (6, False)])
def test_native_layout_needs_no_swaps(backend, n, periodic):
    layout = hw.find_layout(backend, n, periodic)
    qc = tfim_circuit(n, 1.0, 2.0, 0.25, 4, periodic=periodic)
    rep = hw.transpile_report(qc, hw.to_isa(qc, backend, layout))
    assert rep["swaps"] == 0 and rep["extra_two_qubit"] == 0


def test_six_site_ring_is_rejected(backend):
    with pytest.raises(ValueError, match="no cycle of 6"):
        hw.find_layout(backend, 6, periodic=True)


def test_folding_scales_two_qubit_gates_exactly(backend):
    qc = tfim_circuit(6, 1.0, 1.0, 0.25, 2, periodic=False)
    isa = hw.to_isa(qc, backend, hw.find_layout(backend, 6, periodic=False))
    _, factors = hw.zne_circuits(isa, 6, backend, (1, 3, 5))
    assert factors == [1.0, 3.0, 5.0]


def test_ideal_counts_match_statevector(backend):
    n, periodic = 6, False
    qc = tfim_circuit(n, 1.0, 1.0, 0.25, 3, periodic=periodic)
    isa = hw.to_isa(qc, backend, hw.find_layout(backend, n, periodic))
    circs, _ = hw.zne_circuits(isa, n, backend, (1, 3))
    counts = hw.run_simulator(circs, simulator(mode="ideal"), shots=20000)
    exact = state_expectations(trotter_state(n, 1.0, 1.0, 0.75, 3, periodic=periodic), n, periodic)
    for c in counts:  # folding is the identity without noise
        got = counts_expectations(c, n, periodic)
        for k in ("Zavg", "Mzz"):
            assert got[k] == pytest.approx(exact[k], abs=0.02)
