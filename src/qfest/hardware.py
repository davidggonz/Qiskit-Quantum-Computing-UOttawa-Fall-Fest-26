"""IBM backend setup, layout, transpilation, ZNE folding and job submission.  (Owner: Boutaina)

Ported from Boutaina's `hardware.py` (starswithboutaina/Qiskit-Fall-Fest-2026-uOttawa-UO-Qubit,
branch feature/canonical-tfim), adapted to the shared `qfest` modules, with these fixes:

  * Layout: the logical ring/chain is placed on a physical cycle/path of the coupling map
    (`find_layout`), so transpilation adds zero SWAPs. A periodic ring needs a cycle of the
    same length: on heavy-hex the shortest is 12, so N = 6 periodic is rejected (use N = 12,
    or an open chain).
  * ZNE folding is done on the transpiled (ISA) circuit and re-transpiled at optimization
    level 0 on the same physical qubits, so the 2-qubit gate count scales exactly 1:3:5 and
    every noise factor runs on the same qubits (previously each scale got a new layout and
    the effective factors were 1:3.3:5.4).
  * Hardware jobs receive ISA circuits (Runtime rejects untranspiled circuits) and run in job
    or batch mode via SamplerV2; observables are computed from counts with `qfest.shots`, the
    same code path as the simulator.
  * No silent fallbacks: if Aer or Runtime fails, the error is raised.

Credentials: set QISKIT_IBM_TOKEN (and optionally QISKIT_IBM_INSTANCE) in the environment or in
a local `.env` file (git-ignored), or use QiskitRuntimeService.save_account. Never commit tokens.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import rustworkx as rx
from qiskit import ClassicalRegister, transpile

TWO_QUBIT_GATES = ("cz", "ecr", "cx", "rzz")
FAKE_BACKENDS = ("marrakesh", "fez", "torino", "kingston")


# ---------------------------------------------------------------------------
# Backend
# ---------------------------------------------------------------------------

def _load_env_file(path=".env"):
    """Load KEY=VALUE lines from a local .env file into os.environ (existing vars win)."""
    p = Path(path)
    if not p.exists():
        return
    for line in p.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, v = line.split("=", 1)
            os.environ.setdefault(k.strip(), v.strip().strip("'\""))


def get_service(channel="ibm_quantum_platform"):
    from qiskit_ibm_runtime import QiskitRuntimeService

    _load_env_file()
    token = os.environ.get("QISKIT_IBM_TOKEN", "").strip()
    instance = os.environ.get("QISKIT_IBM_INSTANCE", "").strip() or None
    if token:
        return QiskitRuntimeService(channel=channel, token=token, instance=instance)
    return QiskitRuntimeService()  # falls back to a saved account


def get_backend(name=None, fake=None, min_qubits=12):
    """A fake backend (`fake="marrakesh"`, offline) or a real one (by name, else least busy)."""
    if fake:
        from qiskit_ibm_runtime import fake_provider

        cls = getattr(fake_provider, "Fake" + fake.removeprefix("ibm_").capitalize())
        return cls()
    service = get_service()
    if name:
        return service.backend(name)
    return service.least_busy(operational=True, simulator=False, min_num_qubits=min_qubits)


def _graph(backend):
    return backend.coupling_map.graph.to_undirected(multigraph=False)


def _two_qubit_gate(backend):
    names = set(backend.target.operation_names)
    for g in TWO_QUBIT_GATES:
        if g in names:
            return g
    raise ValueError(f"no known 2-qubit gate in {sorted(names)}")


def _edge_error(backend, a, b):
    props = backend.target[_two_qubit_gate(backend)]
    p = props.get((a, b)) or props.get((b, a))
    return p.error if p is not None and p.error is not None else 1.0


def _readout_error(backend, q):
    p = backend.target["measure"].get((q,))
    return p.error if p is not None and p.error is not None else 1.0


def layout_cost(backend, qubits, periodic):
    """Sum of 2-qubit errors on the used edges plus readout errors (lower is better)."""
    n = len(qubits)
    edges = [(qubits[i], qubits[i + 1]) for i in range(n - 1)]
    if periodic:
        edges.append((qubits[-1], qubits[0]))
    return sum(_edge_error(backend, a, b) for a, b in edges) + sum(_readout_error(backend, q) for q in qubits)


def _ordered_cycles(g, n):
    for cyc in rx.cycle_basis(g):
        if len(cyc) == n and all(g.has_edge(cyc[i], cyc[(i + 1) % n]) for i in range(n)):
            yield list(cyc)


def _paths(g, n, max_paths=2000):
    """Enumerate up to `max_paths` simple paths with n nodes (DFS)."""
    found = []

    def dfs(path):
        if len(found) >= max_paths:
            return
        if len(path) == n:
            found.append(list(path))
            return
        for nb in g.neighbors(path[-1]):
            if nb not in path:
                path.append(nb)
                dfs(path)
                path.pop()

    for start in g.node_indices():
        dfs([start])
        if len(found) >= max_paths:
            break
    return found


def find_layout(backend, n, periodic=True):
    """Physical qubits for a ring (cycle) or chain (path) of n qubits, lowest `layout_cost` first.

    Logical qubit i -> returned[i]; consecutive entries (and last->first for a ring) are coupled.
    """
    g = _graph(backend)
    cands = list(_ordered_cycles(g, n)) if periodic else _paths(g, n)
    if not cands:
        lengths = sorted({len(c) for c in rx.cycle_basis(g)})
        raise ValueError(
            f"no {'cycle' if periodic else 'path'} of {n} qubits on {backend.name}"
            + (f" (cycle lengths available: {lengths[:5]}); use n in those, or an open chain" if periodic else "")
        )
    return min(cands, key=lambda q: layout_cost(backend, q, periodic))


def embed_ring(n, backend):
    """Brief-compatible alias: best physical cycle for a periodic ring of n qubits."""
    return find_layout(backend, n, periodic=True)


# ---------------------------------------------------------------------------
# Transpilation
# ---------------------------------------------------------------------------

def to_isa(qc, backend, layout, optimization_level=1, seed=42):
    """Transpile a measurement-free logical circuit onto `layout` (no routing expected)."""
    if qc.num_clbits:
        raise ValueError("pass the circuit without measurements; use add_measurements()")
    return transpile(qc, backend=backend, initial_layout=layout,
                     optimization_level=optimization_level, seed_transpiler=seed)


def physical_qubits(isa):
    """Physical qubit holding logical qubit i at the end of the circuit."""
    return list(isa.layout.final_index_layout())


def add_measurements(circ, qubits):
    """Measure physical qubits[i] (= logical qubit i) into classical bit i."""
    out = circ.copy()
    creg = ClassicalRegister(len(qubits), "meas")
    out.add_register(creg)
    for i, q in enumerate(qubits):
        out.measure(q, creg[i])
    return out


def two_qubit_count(circ):
    ops = circ.count_ops()
    return int(sum(ops.get(g, 0) for g in TWO_QUBIT_GATES))


def transpile_report(logical, isa):
    """Depth / gate counts of the ISA circuit and SWAP overhead vs the logical circuit.

    With a native layout each logical RZZ costs 2 CZ/ECR, so extra = 2q - 2 * n_rzz should be 0.
    """
    n_rzz = int(logical.count_ops().get("rzz", 0))
    twoq = two_qubit_count(isa)
    ops = isa.count_ops()
    return {
        "depth": int(isa.depth()),
        "size": int(isa.size()),
        "two_qubit": twoq,
        "logical_rzz": n_rzz,
        "extra_two_qubit": twoq - 2 * n_rzz,
        "swaps": int(ops.get("swap", 0)),
        "counts": {str(k): int(v) for k, v in ops.items()},
    }


def transpile_comparison(logical, backend, layout, levels=(0, 1, 2, 3), seed=42):
    """Boutaina's optimization-level comparison, now on a fixed layout and a fixed seed."""
    return [{"optimization_level": lvl, **transpile_report(logical, to_isa(logical, backend, layout, lvl, seed))}
            for lvl in levels]


# ---------------------------------------------------------------------------
# ZNE folding on the ISA circuit
# ---------------------------------------------------------------------------

def fold_isa(isa, scale, backend):
    """Global folding U (U^dag U)^((scale-1)/2) on the physical circuit.

    Re-transpiled at optimization level 0 with a trivial layout so the folds are neither
    cancelled nor moved to other qubits.
    """
    if scale < 1 or scale % 2 == 0:
        raise ValueError("scale must be an odd integer >= 1")
    folded = isa.copy()
    inv = isa.inverse()
    for _ in range((scale - 1) // 2):
        folded = folded.compose(inv).compose(isa)
    return transpile(folded, backend=backend, optimization_level=0,
                     initial_layout=list(range(folded.num_qubits)))


def zne_circuits(isa, n, backend, scales=(1, 3, 5)):
    """Measured folded circuits and their effective noise factors (2-qubit gate ratios)."""
    base = two_qubit_count(isa)
    qubits = physical_qubits(isa)[:n]
    circs, factors = [], []
    for s in scales:
        f = fold_isa(isa, s, backend)
        circs.append(add_measurements(f, qubits))
        factors.append(two_qubit_count(f) / base if base else float(s))
    return circs, factors


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

def run_simulator(circuits, simulator, shots=4000, seed=42):
    """Counts from an AerSimulator (see qfest.noise.simulator). Errors are raised, not hidden."""
    res = simulator.run(list(circuits), shots=shots, seed_simulator=seed).result()
    return [dict(res.get_counts(i)) for i in range(len(circuits))]


def run_batch(circuits, backend, shots=4000, dynamical_decoupling=True, twirling=True, batch=False):
    """Run ISA circuits on IBM hardware with SamplerV2 (job mode, or one Batch if batch=True).

    Returns (counts list, job id). Dynamical decoupling and gate twirling are Sampler options;
    readout (TREX) mitigation needs the Estimator and is not applied here.
    """
    from qiskit_ibm_runtime import Batch, SamplerV2

    def _run(mode):
        sampler = SamplerV2(mode=mode)
        sampler.options.default_shots = shots
        sampler.options.dynamical_decoupling.enable = dynamical_decoupling
        sampler.options.twirling.enable_gates = twirling
        job = sampler.run(list(circuits))
        result = job.result()
        return [r.data.meas.get_counts() for r in result], job.job_id()

    if batch:
        with Batch(backend=backend) as b:
            return _run(b)
    return _run(backend)
