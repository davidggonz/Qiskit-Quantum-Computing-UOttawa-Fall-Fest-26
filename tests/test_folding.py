import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # hardware.py lives at the repo root

from qiskit import transpile  # noqa: E402
from qiskit.quantum_info import Operator  # noqa: E402

import hardware as hw  # noqa: E402


@pytest.fixture(scope="module")
def base():
    qc = hw.build_tfim_circuit(4, 1.0, 1.0, 0.1, 3, add_measure=False)
    return transpile(qc, basis_gates=hw.SIM_BASIS, optimization_level=1, seed_transpiler=42)


@pytest.mark.parametrize("scale", [1, 2, 3, 4, 5])
def test_folding_preserves_unitary(base, scale):
    folded = hw.fold_two_qubit_gates(base, scale)
    assert Operator(folded).equiv(Operator(base))


@pytest.mark.parametrize("scale", [1, 2, 3, 4, 5])
def test_folding_scales_two_qubit_count(base, scale):
    folded = hw.fold_two_qubit_gates(base, scale)
    n2q = sum(1 for i in base.data if len(i.qubits) == 2)
    # partial folds add one pair per selected gate, so the scale is exact to within 2/n2q
    assert hw.effective_scale(base, folded) == pytest.approx(scale, abs=2.0 / n2q)
