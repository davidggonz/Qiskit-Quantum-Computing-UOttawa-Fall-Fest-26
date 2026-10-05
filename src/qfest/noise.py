"""Local noise models for TFIM simulation.

The upstream qfest-plots source was unavailable. ``simple_model`` therefore
uses the project's pre-existing 1%/3% depolarizing-noise demo parameters.
"""

from qiskit_aer.noise import NoiseModel, depolarizing_error


def simple_model() -> NoiseModel:
    """Return a nontrivial depolarizing model for common IBM basis gates."""
    model = NoiseModel()
    model.add_all_qubit_quantum_error(
        depolarizing_error(0.01, 1), ["id", "rz", "sx", "x", "h", "rx"]
    )
    model.add_all_qubit_quantum_error(
        depolarizing_error(0.03, 2), ["cx", "cz", "ecr", "rzz"]
    )
    return model


def from_backend(backend) -> NoiseModel:
    """Build an Aer noise model from a backend's calibrated properties."""
    from qiskit_aer import AerSimulator

    model = AerSimulator.from_backend(backend).options.noise_model
    if model is None or not model.to_dict().get("errors"):
        raise ValueError("backend did not provide a nontrivial noise model")
    return model
