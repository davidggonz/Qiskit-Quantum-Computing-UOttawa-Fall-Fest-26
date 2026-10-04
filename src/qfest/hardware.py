"""IBM backend setup, transpilation and job submission.  (Owner: Boutaina)

TODO:
  - get_backend(name=None): QiskitRuntimeService, least-busy Heron/Eagle, or a fake backend
    (qiskit_ibm_runtime.fake_provider) when offline.
  - embed_ring(n, backend): choose a physical qubit cycle of length n in the coupling map
    (initial_layout) so a periodic ring needs zero SWAPs.
  - transpile_report(qc, backend): return 2q gate count, depth, SWAP count.
  - run_batch(circuits, observables, backend, mitigation): Estimator with twirling, DD, TREX, ZNE
    options; submit as ONE batch/session to save QPU time; save raw job ids to results/.
Never commit API tokens. Use QiskitRuntimeService.save_account locally or an env var.
"""


def get_backend(name=None):
    raise NotImplementedError


def embed_ring(n, backend):
    raise NotImplementedError


def transpile_report(qc, backend):
    raise NotImplementedError


def run_batch(circuits, observables, backend, mitigation):
    raise NotImplementedError
