"""Hardware-side module for TFIM Trotter simulation on IBM Quantum hardware.

Covers: backend setup, TFIM Trotter circuits, transpilation comparison,
hardware-aware layout selection, ZNE error mitigation, Aer + hardware
execution, exact-diagonalization analysis, JSON + PNG outputs.

Physics follows the team contract (TEAM.md) and the previous report:
  - circuit: `qfest.tfim.tfim_circuit` (2nd-order Trotter, fused Rx, periodic ring);
  - observables: Mz = sqrt(<(sum Z)^2>)/N (RMS magnitude, NOT the signed mean),
    Mx = <sum X>/N, Mzz = <sum_bonds Z Z>/N, bonds matching the circuit;
  - exact reference: `qfest.ed.quench` (exact evolution, same boundary conditions).

Compatible with Qiskit 1.0+ (tested on Qiskit 2.5.2,
qiskit-aer 0.17.2, qiskit-ibm-runtime 0.50.0).

Secrets: never hardcode tokens. Copy `.env.example` to `.env`:
    QISKIT_IBM_TOKEN=<PINQ2 or personal token>
    QISKIT_IBM_CHANNEL=ibm_cloud
    QISKIT_IBM_INSTANCE=<CRN or empty for auto>
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

# Use the team package (src/qfest) even without `pip install -e .` (e.g. Colab upload).
_SRC = Path(__file__).resolve().parent / "src"
if _SRC.exists() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from qfest import ed as qed  # noqa: E402
from qfest.extrapolation import EXTRAPOLATORS  # noqa: E402
from qfest.tfim import tfim_circuit  # noqa: E402

OBS_KEYS = ("Mz", "Mx", "Mzz")

# ---------------------------------------------------------------------------
# .env loader (no extra dependency)
# ---------------------------------------------------------------------------

def load_dotenv(path: str | Path = ".env") -> None:
    """Load KEY=VALUE lines from .env into os.environ (no override)."""
    p = Path(path)
    if not p.exists():
        return
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip().strip('"').strip("'")
        os.environ.setdefault(k, v)


load_dotenv()

DEFAULT_CHANNEL = os.getenv("QISKIT_IBM_CHANNEL", "ibm_cloud")
PREFERRED_BACKENDS = ["ibm_marrakesh", "ibm_fez", "ibm_torino", "ibm_kyiv"]


# ---------------------------------------------------------------------------
# 1. Backend setup
# ---------------------------------------------------------------------------

def get_service(token: Optional[str] = None, channel: Optional[str] = None,
                instance: Optional[str] = None):
    """Connect to IBM Quantum via QiskitRuntimeService.

    Returns None in offline mode (no token) so the rest still runs on Aer.
    """
    token = token or os.getenv("QISKIT_IBM_TOKEN", "")
    channel = channel or os.getenv("QISKIT_IBM_CHANNEL", DEFAULT_CHANNEL)
    instance = instance or os.getenv("QISKIT_IBM_INSTANCE", "") or None
    if not token:
        print("[backend] No QISKIT_IBM_TOKEN found — offline/Aer mode.")
        return None
    from qiskit_ibm_runtime import QiskitRuntimeService
    kwargs: Dict[str, Any] = {"channel": channel, "token": token}
    if instance:
        kwargs["instance"] = instance
    service = QiskitRuntimeService(**kwargs)
    print(f"[backend] Connected (channel={channel}).")
    return service


def select_backend(service, preferred: Sequence[str] = tuple(PREFERRED_BACKENDS)):
    """Select a real backend, preferring `preferred` names, else least-busy."""
    backends = service.backends(operational=True, simulator=False)
    by_name = {b.name: b for b in backends}
    for name in preferred:
        if name in by_name:
            print(f"[backend] Selected preferred backend: {name}")
            return by_name[name]
    cand = service.least_busy(operational=True, simulator=False)
    print(f"[backend] Preferred not found, using least busy: {cand.name}")
    return cand


def get_fake_backend(name: str = "marrakesh"):
    """Offline stand-in with a real Heron coupling map and calibration (no QPU time)."""
    from qiskit_ibm_runtime import fake_provider
    cls = getattr(fake_provider, f"Fake{name.capitalize()}")
    return cls()


def backend_properties_dict(backend) -> Dict[str, Any]:
    """Collect printable backend properties + per-qubit/gate error rates."""
    target = getattr(backend, "target", None)
    n = backend.num_qubits
    try:
        cmap = sorted([list(e) for e in backend.coupling_map.get_edges()]) \
            if backend.coupling_map else []
    except Exception:
        cmap = []
    try:
        basis = sorted(list(target.operation_names)) if target else []
    except Exception:
        basis = []

    readout_err: Dict[int, float] = {}
    err_1q: Dict[int, float] = {}
    err_2q: Dict[Tuple[int, int], float] = {}
    if target:
        for q in range(n):
            # Physical 1q gates only: rz is virtual (error 0) and id is a delay,
            # so including them made every median 1q error read 0.0.
            errs = []
            for op in ("x", "sx"):
                try:
                    e = target[op][(q,)].error
                    if e is not None:
                        errs.append(float(e))
                except Exception:
                    continue
            if errs:
                err_1q[q] = float(np.mean(errs))
            try:
                readout_err[q] = float(target["measure"][(q,)].error)
            except Exception:
                pass
        for op in ("cx", "cz", "ecr"):
            try:
                qargs = target.qargs_for_operation_name(op)
            except Exception:
                continue
            for qa in qargs:
                if len(qa) == 2:
                    try:
                        err_2q[tuple(qa)] = float(target[op][qa].error)
                    except Exception:
                        pass

    info = {
        "name": getattr(backend, "name", str(backend)),
        "num_qubits": n,
        "coupling_map": cmap,
        "basis_gates": basis[:24],
        "dt": getattr(backend, "dt", None),
        "readout_err": {str(k): v for k, v in readout_err.items()},
        "median_1q_err": float(np.median(list(err_1q.values()))) if err_1q else None,
        "median_2q_err": float(np.median(list(err_2q.values()))) if err_2q else None,
        "median_readout_err": float(np.median(list(readout_err.values()))) if readout_err else None,
    }
    print(f"[backend] {info['name']}: {n} qubits, "
          f"median 1q={info['median_1q_err']}, "
          f"2q={info['median_2q_err']}, ro={info['median_readout_err']}")
    print(f"[backend] basis (subset): {info['basis_gates']}")
    print(f"[backend] coupling edges: {len(cmap)}")
    return info


# ---------------------------------------------------------------------------
# 2. TFIM Trotter circuit and observables
# ---------------------------------------------------------------------------

def build_tfim_circuit(n: int = 4, j: float = 1.0, h: float = 0.5,
                       dt: float = 0.1, steps: int = 3,
                       periodic: bool = True, order: int = 2,
                       add_measure: bool = True, basis: str = "z"):
    """Trotterized TFIM, H = -J Σ ZZ - h Σ X, quench from |0>^N.

    Delegates to the team's `qfest.tfim.tfim_circuit` (2nd order by default,
    fused Rx half-steps, edge-coloured ZZ layers). `basis="x"` adds H gates
    before measuring, needed for Mx from counts.

    order=1 (ZZ layer then X layer) is NOT worse for Mz/Mzz from |0...0>: that
    state is a ZZ eigenstate, so the product equals a ZZ-outer symmetric formula
    up to a diagonal phase. It is clearly worse for Mx (see
    experiments/compare_hardware_fix.py, figure 2).
    """
    qc = tfim_circuit(n, j, h, dt, steps, order=order, periodic=periodic)
    if add_measure:
        if basis == "x":
            qc.h(range(n))
        qc.measure_all()
    qc.metadata = {"n": n, "J": j, "h": h, "dt": dt, "steps": steps,
                   "periodic": periodic, "order": order, "basis": basis}
    return qc


def _pauli(n: int, ops: Dict[int, str]) -> str:
    """Qiskit label (qubit 0 is the rightmost character)."""
    label = ["I"] * n
    for q, p in ops.items():
        label[n - 1 - q] = p
    return "".join(label)


def build_tfim_observables(n: int = 4, j: float = 1.0, h: float = 0.5,
                           periodic: bool = True):
    """SparsePauliOps for the energy and the report's observables.

    Bonds come from `qfest.ed.bonds`, so the operators always match the circuit.
    Mz itself is not linear in the state: we return Mz2 = (Σ Z)^2 / N^2 and the
    caller takes Mz = sqrt(<Mz2>) (see `observables_from_evs`).
    """
    from qiskit.quantum_info import SparsePauliOp
    bnds = qed.bonds(n, periodic)
    energy = SparsePauliOp([_pauli(n, {a: "Z", b: "Z"}) for a, b in bnds], [-j] * len(bnds)) \
        + SparsePauliOp([_pauli(n, {q: "X"}) for q in range(n)], [-h] * n)
    mx = SparsePauliOp([_pauli(n, {q: "X"}) for q in range(n)], [1.0 / n] * n)
    mzz = SparsePauliOp([_pauli(n, {a: "Z", b: "Z"}) for a, b in bnds], [1.0 / n] * len(bnds))
    # (Σ Z)^2 = N + 2 Σ_{a<b} Z_a Z_b
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    mz2 = SparsePauliOp(["I" * n] + [_pauli(n, {a: "Z", b: "Z"}) for a, b in pairs],
                        [n / n ** 2] + [2.0 / n ** 2] * len(pairs))
    return {"energy": energy.simplify(), "Mz2": mz2.simplify(), "Mx": mx, "Mzz": mzz}


def observables_from_evs(evs: Dict[str, float]) -> Dict[str, float]:
    """Map raw expectation values {Mz2, Mx, Mzz} to the report's {Mz, Mx, Mzz}."""
    return {"Mz": float(np.sqrt(max(evs["Mz2"], 0.0))), "Mx": float(evs["Mx"]),
            "Mzz": float(evs["Mzz"])}


def exact_tfim_observables(n: int = 4, j: float = 1.0, h: float = 0.5,
                           t: float = 0.3, periodic: bool = True) -> Dict[str, float]:
    """Exact (not Trotterized) Mz, Mx, Mzz at time t via `qfest.ed.quench`."""
    r = qed.quench(n, [t], j, h, periodic)
    return {k: float(r[k][0]) for k in OBS_KEYS}


def exact_tfim_magnetization(n: int = 4, j: float = 1.0, h: float = 0.5,
                             dt: float = 0.1, steps: int = 3,
                             periodic: bool = True) -> float:
    """Exact Mz (RMS magnitude) at t = dt*steps. Kept for the old call signature."""
    return exact_tfim_observables(n, j, h, dt * steps, periodic)["Mz"]


# ---------------------------------------------------------------------------
# 3. Transpilation comparison
# ---------------------------------------------------------------------------

def estimate_circuit_error(circuit, backend) -> Optional[float]:
    """1 - Π(1-e_g) using backend.target error rates; None if unavailable."""
    target = getattr(backend, "target", None)
    if target is None:
        return None
    qc = circuit.remove_final_measurements(inplace=False)
    fid = 1.0
    for inst in qc.data:
        qa = tuple(qc.find_bit(q).index for q in inst.qubits)
        try:
            e = target[inst.operation.name][qa].error
            if e is not None:
                fid *= (1.0 - float(e))
        except Exception:
            continue
    return float(1.0 - fid)


def two_qubit_count(counts: Dict[str, int]) -> int:
    return int(sum(v for k, v in counts.items() if k in ("cx", "cz", "ecr", "rzz", "swap")))


def transpile_compare(circuit, backend=None, levels: Sequence[int] = (0, 1, 2, 3),
                      initial_layout: Optional[List[int]] = None):
    """Transpile at each opt level; record depth, gates, 2q count, est. error."""
    from qiskit import transpile
    rows = []
    for lv in levels:
        t = transpile(circuit, backend=backend, optimization_level=lv,
                      initial_layout=initial_layout, seed_transpiler=42)
        counts = dict(t.count_ops())
        row = {
            "optimization_level": lv,
            "depth": int(t.depth()),
            "qubits": int(t.num_qubits),
            "total_gates": int(sum(counts.values())),
            "cx_like": two_qubit_count(counts),
            "counts": {k: int(v) for k, v in counts.items()},
            "est_error": estimate_circuit_error(t, backend) if backend else None,
        }
        rows.append(row)
    print(f"{'level':>5} {'depth':>6} {'gates':>6} {'2q':>5} {'est_err':>9}")
    for r in rows:
        e = f"{r['est_error']:.4f}" if r["est_error"] is not None else "n/a"
        print(f"{r['optimization_level']:>5} {r['depth']:>6} "
              f"{r['total_gates']:>6} {r['cx_like']:>5} {e:>9}")
    return rows


# ---------------------------------------------------------------------------
# 4. Hardware-aware layout
# ---------------------------------------------------------------------------

def _error_tables(backend):
    target = backend.target
    q_score = {}
    for q in range(backend.num_qubits):
        s = 0.0
        try:
            s += float(target["measure"][(q,)].error)
        except Exception:
            s += 0.02
        errs = []
        for op in ("x", "sx"):
            try:
                errs.append(float(target[op][(q,)].error))
            except Exception:
                pass
        s += float(np.mean(errs)) if errs else 0.001
        q_score[q] = s

    def edge_err(a, b):
        best = None
        for op in ("cx", "cz", "ecr"):
            for qa in ((a, b), (b, a)):
                try:
                    e = float(target[op][qa].error)
                    best = e if best is None else min(best, e)
                except Exception:
                    pass
        return best if best is not None else 0.02

    neighbors = {q: set(backend.coupling_map.neighbors(q)) for q in range(backend.num_qubits)}
    return q_score, edge_err, neighbors


def _path_score(path, q_score, edge_err, closed):
    s = sum(q_score[q] for q in path) + sum(edge_err(a, b) for a, b in zip(path, path[1:]))
    return s + (edge_err(path[-1], path[0]) if closed else 0.0)


def find_rings(backend, n: int, limit: int = 2000) -> List[List[int]]:
    """All simple cycles of exactly n physical qubits (each found once).

    Heavy-hex has no cycles shorter than 12, so N = 6/8/10 rings cannot be
    embedded without SWAPs; N = 12 (one hexagon) can.
    """
    _, _, nb = _error_tables(backend)
    rings = []

    def dfs(start, path, seen):
        if len(rings) >= limit:
            return
        last = path[-1]
        if len(path) == n:
            if start in nb[last] and path[1] < path[-1]:  # drop the reversed duplicate
                rings.append(list(path))
            return
        for q in nb[last]:
            if q > start and q not in seen:
                seen.add(q)
                path.append(q)
                dfs(start, path, seen)
                path.pop()
                seen.discard(q)

    for s in range(backend.num_qubits):
        dfs(s, [s], {s})
    return rings


def find_chain(backend, n: int) -> List[int]:
    """Lowest-error simple path of n physical qubits, in path order
    (layout[i] and layout[i+1] are always coupled)."""
    q_score, edge_err, nb = _error_tables(backend)
    best, best_s = None, float("inf")
    for start in sorted(q_score, key=q_score.get)[:40]:
        path = [start]
        while len(path) < n:
            cand = [q for q in nb[path[-1]] if q not in path]
            if not cand:
                break
            path.append(min(cand, key=lambda q: edge_err(path[-1], q) + q_score[q]))
        if len(path) == n:
            s = _path_score(path, q_score, edge_err, False)
            if s < best_s:
                best, best_s = path, s
    if best is None:
        raise RuntimeError(f"no simple path of {n} qubits found")
    return best


def select_best_qubits(backend, n: int, periodic: bool = True) -> List[int]:
    """Physical qubits for the TFIM, ordered so the circuit's bonds are couplers.

    periodic: lowest-error n-cycle (0 SWAPs). Falls back to a chain, with a
    warning, when no n-cycle exists (the closing bond then costs SWAPs).
    """
    q_score, edge_err, _ = _error_tables(backend)
    if periodic:
        rings = find_rings(backend, n)
        if rings:
            best = min(rings, key=lambda r: _path_score(r, q_score, edge_err, True))
            print(f"[layout] best {n}-ring of {len(rings)}: {best}")
            return best
        print(f"[layout] WARNING: no {n}-qubit cycle on {backend.name}; "
              f"using a chain, the closing bond will need SWAPs")
    chain = find_chain(backend, n)
    print(f"[layout] best {n}-chain: {chain}")
    return chain


def hardware_aware_transpile(circuit, backend, n: int, periodic: bool = True):
    """Transpile opt_level=3 with initial_layout=best qubits; compare to default."""
    from qiskit import transpile
    layout = select_best_qubits(backend, n, periodic)
    default = transpile(circuit, backend=backend, optimization_level=1, seed_transpiler=42)
    aware = transpile(circuit, backend=backend, optimization_level=3,
                      initial_layout=layout, seed_transpiler=42)
    for label, t in (("default(lv1)", default), ("aware(lv3+layout)", aware)):
        c = t.count_ops()
        print(f"[layout] {label}: depth={t.depth()}, 2q={two_qubit_count(c)}, {dict(c)}")
    return aware, {"layout": layout,
                   "default": {"depth": default.depth(), "counts": dict(default.count_ops())},
                   "aware": {"depth": aware.depth(), "counts": dict(aware.count_ops())}}


# ---------------------------------------------------------------------------
# 5-6. Execution (Aer + hardware) and ZNE
# ---------------------------------------------------------------------------

SIM_BASIS = ["cz", "rz", "sx", "x"]  # Heron-like basis for Aer runs


def noise_model_like(backend=None, p1q: float = 3e-4, p2q: float = 3e-3,
                     p_meas: float = 1e-2):
    """Depolarizing + readout model on SIM_BASIS; medians from `backend` if given.

    Lightweight stand-in until Toto's `qfest.noise` lands. Uses backend medians
    rather than NoiseModel.from_backend so it stays small (no 156-qubit layout).
    """
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error
    if backend is not None:
        info = backend_properties_dict(backend)
        p1q = info["median_1q_err"] or p1q
        p2q = info["median_2q_err"] or p2q
        p_meas = info["median_readout_err"] or p_meas
    nm = NoiseModel(basis_gates=SIM_BASIS)
    nm.add_all_qubit_quantum_error(depolarizing_error(p1q, 1), ["sx", "x"])
    nm.add_all_qubit_quantum_error(depolarizing_error(p2q, 2), ["cz"])
    nm.add_all_qubit_readout_error(ReadoutError([[1 - p_meas, p_meas], [p_meas, 1 - p_meas]]))
    return nm


def run_counts_aer(circuit, shots: int = 4096, seed: int = 42,
                   noise_model=None) -> Dict[str, int]:
    """Counts on AerSimulator (ideal unless `noise_model` is given)."""
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    sim = AerSimulator(seed_simulator=seed, noise_model=noise_model)
    qc = transpile(circuit, basis_gates=SIM_BASIS, optimization_level=1, seed_transpiler=42) \
        if noise_model is not None else circuit
    job = sim.run(qc, shots=shots)
    return {k: int(v) for k, v in job.result().get_counts().items()}


def run_expectation_aer(circuit, observables: Dict[str, Any], shots: int = 4096,
                        seed: int = 42, noise_model=None,
                        optimization_level: int = 1) -> Dict[str, float]:
    """Raw expectation values {name: <O>} via Aer EstimatorV2.

    Pass optimization_level=0 for folded circuits, otherwise the transpiler
    cancels the G·G†·G folds and ZNE sees no extra noise.
    """
    from qiskit import transpile
    from qiskit_aer.primitives import EstimatorV2 as AerEstimator
    qc = circuit.remove_final_measurements(inplace=False)
    opts: Dict[str, Any] = {"default_precision": 1.0 / np.sqrt(shots),
                            "run_options": {"seed": seed}}
    if noise_model is not None:
        opts["backend_options"] = {"noise_model": noise_model}
    est = AerEstimator(options=opts)
    tqc = transpile(qc, basis_gates=SIM_BASIS, optimization_level=optimization_level,
                    seed_transpiler=42)
    names = list(observables)
    res = est.run([(tqc, [observables[k] for k in names])]).result()
    return {k: float(np.real(v)) for k, v in zip(names, res[0].data.evs)}


def zne_extrapolate(noise_factors: Sequence[float], values: Sequence[float],
                    method: str = "linear") -> float:
    """Extrapolate to zero noise with the team's `qfest.extrapolation`
    (linear | richardson | exponential; exponential fits a*exp(-b x) + c)."""
    return float(EXTRAPOLATORS[method](list(noise_factors), list(values)))


def fold_two_qubit_gates(circuit, scale: int):
    """Local unitary folding: every 2q gate G -> G (G† G)^((scale-1)/2).

    Apply AFTER transpiling to the target basis (cz/ecr/cx/rzz); before
    transpiling the TFIM circuit only has rzz, and repeating a non-self-inverse
    rzz would change the unitary. scale=1 -> identity.
    """
    if scale == 1:
        return circuit.copy()
    assert scale % 2 == 1, "scale must be odd"
    folded = circuit.copy_empty_like()
    for inst in circuit.data:
        folded.append(inst)
        if len(inst.qubits) == 2 and inst.operation.name in ("cx", "cz", "ecr", "rzz"):
            for _ in range((scale - 1) // 2):
                folded.append(inst.operation.inverse(), inst.qubits)
                folded.append(inst.operation, inst.qubits)
    folded.metadata = dict(getattr(circuit, "metadata", {}) or {})
    return folded


def run_with_zne_aer(circuit, observables: Dict[str, Any], shots: int = 4096,
                     seed: int = 42, noise_model=None,
                     noise_factors: Sequence[int] = (1, 3, 5),
                     extrapolators: Sequence[str] = ("linear", "richardson", "exponential"),
                     ) -> Dict[str, Any]:
    """Manual ZNE on Aer: transpile, fold 2q gates, evaluate, extrapolate.

    Needs a noise model: on the ideal simulator every factor gives the same
    value and ZNE is a no-op.
    """
    from qiskit import transpile
    if noise_model is None:
        print("[zne-aer] WARNING: no noise model, all noise factors give the same value")
    base = transpile(circuit.remove_final_measurements(inplace=False),
                     basis_gates=SIM_BASIS, optimization_level=1, seed_transpiler=42)
    raw = {k: [] for k in OBS_KEYS}
    for f in noise_factors:
        evs = run_expectation_aer(fold_two_qubit_gates(base, int(f)), observables,
                                  shots=shots, seed=seed, noise_model=noise_model,
                                  optimization_level=0)
        for k, v in observables_from_evs(evs).items():
            raw[k].append(v)
    mitigated = {m: {k: zne_extrapolate(noise_factors, raw[k], m) for k in OBS_KEYS}
                 for m in extrapolators}
    print(f"[zne-aer] raw={raw}")
    print(f"[zne-aer] mitigated={mitigated}")
    return {"noise_factors": list(noise_factors), "raw": raw, "mitigated": mitigated}


def run_expectation_hardware(circuit, observables: Dict[str, Any], backend,
                             shots: int = 4096, resilience_level: int = 1,
                             use_zne: bool = False, initial_layout=None) -> Dict[str, Any]:
    """Expectation values on a real backend: ONE EstimatorV2 job, all observables.

    Job mode (mode=backend), no Session: Sessions are not available on the
    Open plan, and one job with one PUB is the cheapest submission.
    """
    from qiskit import transpile
    from qiskit_ibm_runtime import EstimatorV2
    qc = circuit.remove_final_measurements(inplace=False)
    qc_t = transpile(qc, backend=backend, optimization_level=3,
                     initial_layout=initial_layout, seed_transpiler=42)
    names = list(observables)
    # Observables are defined on n logical qubits; map them onto the ISA layout.
    isa_obs = [observables[k].apply_layout(qc_t.layout) for k in names]
    opts: Dict[str, Any] = {"default_shots": shots, "resilience_level": resilience_level}
    if use_zne:
        opts["resilience_level"] = max(resilience_level, 2)
        opts["resilience"] = {"zne_mitigation": True,
                              "zne": {"noise_factors": (1, 3, 5),
                                      "extrapolator": ("exponential", "linear")}}
    est = EstimatorV2(mode=backend, options=opts)
    job = est.run([(qc_t, isa_obs)])
    res = job.result()
    evs = {k: float(np.real(v)) for k, v in zip(names, res[0].data.evs)}
    out = observables_from_evs(evs)
    print(f"[hw] backend={backend.name} job={job.job_id()} res_lv={opts['resilience_level']} "
          f"zne={use_zne} → {out}")
    return {"values": out, "raw_evs": evs, "job_id": job.job_id(), "shots": shots,
            "options": opts}


def run_counts_hardware(circuit, backend, shots: int = 4096,
                        initial_layout=None) -> Dict[str, int]:
    """Counts on real backend via SamplerV2 (job mode)."""
    from qiskit_ibm_runtime import SamplerV2
    from qiskit import transpile
    qc_t = transpile(circuit, backend=backend, optimization_level=3,
                     initial_layout=initial_layout, seed_transpiler=42)
    samp = SamplerV2(mode=backend, options={"default_shots": shots})
    res = samp.run([(qc_t,)]).result()
    creg = qc_t.cregs[0].name
    counts = getattr(res[0].data, creg).get_counts()
    return {k: int(v) for k, v in counts.items()}


# ---------------------------------------------------------------------------
# 7. Analysis + plots + persistence
# ---------------------------------------------------------------------------

def observables_from_counts(counts: Dict[str, int], n: int,
                            periodic: bool = True) -> Dict[str, float]:
    """Mz (RMS), Mzz and the signed <Z> from Z-basis counts.

    Mz = sqrt(<(Σ z)^2>)/N is averaged per shot BEFORE the square root; the
    signed mean <Σ z>/N is reported separately as "Z_signed" and is NOT Mz.
    """
    total = sum(counts.values())
    bnds = qed.bonds(n, periodic)
    s2 = zz = zs = 0.0
    for bits, c in counts.items():
        b = bits.replace(" ", "")[-n:][::-1]  # b[i] = qubit i
        z = [1.0 if ch == "0" else -1.0 for ch in b]
        s = sum(z)
        s2 += s * s * c
        zs += s * c
        zz += sum(z[a] * z[q] for a, q in bnds) * c
    if not total:
        return {"Mz": 0.0, "Mzz": 0.0, "Z_signed": 0.0}
    return {"Mz": float(np.sqrt(s2 / total) / n), "Mzz": zz / total / n,
            "Z_signed": zs / total / n}


def mx_from_counts(counts: Dict[str, int], n: int) -> float:
    """Mx from counts of the basis="x" circuit."""
    return observables_from_counts(counts, n)["Z_signed"]


def analyze(level_rows, exact: Dict[str, float], aer_vals: Dict[str, Dict[str, float]],
            zne: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Build accuracy table vs the exact baseline, per observable."""
    table = []
    for r in level_rows:
        lv = str(r["optimization_level"])
        v = aer_vals.get(lv)
        entry = {**r, "aer_values": v,
                 "abs_error": {k: abs(v[k] - exact[k]) for k in OBS_KEYS} if v else None}
        table.append(entry)
    return {"exact": exact, "table": table, "zne": zne}


def save_json(payload: Dict[str, Any], path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"[out] JSON → {p}")
    return p


def make_plots(analysis: Dict[str, Any], outdir: str | Path) -> List[Path]:
    """accuracy vs level, depth vs level, error vs mitigation → PNGs."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = analysis["table"]
    levels = [r["optimization_level"] for r in rows]
    colors = {"Mz": "#2a78d6", "Mx": "#eb6834", "Mzz": "#1baf7a"}
    paths = []
    if rows and rows[0].get("abs_error"):
        plt.figure()
        for k in OBS_KEYS:
            plt.plot(levels, [r["abs_error"][k] for r in rows], marker="o", lw=2,
                     color=colors[k], label=k)
        plt.xlabel("optimization level")
        plt.ylabel("|Aer - exact|")
        plt.title("Accuracy vs optimization level")
        plt.legend()
        p = outdir / "accuracy_vs_level.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(p)
    if rows:
        plt.figure()
        plt.plot(levels, [r["depth"] for r in rows], marker="o", lw=2, color="#2a78d6")
        plt.xlabel("optimization level")
        plt.ylabel("transpiled depth")
        plt.title("Depth vs optimization level")
        p = outdir / "depth_vs_level.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(p)
    zne = analysis.get("zne")
    if zne:
        exact = analysis["exact"]
        labels = ["raw (λ=1)"] + list(zne["mitigated"])
        x = np.arange(len(labels))
        plt.figure()
        for i, k in enumerate(OBS_KEYS):
            errs = [abs(zne["raw"][k][0] - exact[k])] + \
                   [abs(zne["mitigated"][m][k] - exact[k]) for m in zne["mitigated"]]
            plt.bar(x + (i - 1) * 0.27, errs, width=0.25, color=colors[k], label=k)
        plt.xticks(x, labels)
        plt.ylabel("absolute error vs exact")
        plt.title("Error: raw vs ZNE extrapolators")
        plt.legend()
        p = outdir / "error_vs_mitigation.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(p)
    for p in paths:
        print(f"[out] PNG → {p}")
    return paths


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="TFIM hardware module driver")
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--J", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.5)
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--order", type=int, default=2, choices=(1, 2))
    ap.add_argument("--open", action="store_true", help="open chain (default: periodic ring)")
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--backend", default="")
    ap.add_argument("--fake", default="marrakesh",
                    help="fake backend for offline transpilation/noise ('' to disable)")
    ap.add_argument("--no-hardware", action="store_true",
                    help="skip real-backend execution (Aer only)")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args(argv)

    seed = 42
    np.random.seed(seed)
    outdir = Path(args.outdir)
    t0 = time.time()
    periodic = not args.open
    t_final = args.dt * args.steps

    qc = build_tfim_circuit(args.n, args.J, args.h, args.dt, args.steps,
                            periodic=periodic, order=args.order)
    core = qc.remove_final_measurements(inplace=False)
    print(f"[circuit] TFIM N={args.n} periodic={periodic} order={args.order} "
          f"depth={core.depth()} gates={dict(core.count_ops())}")
    obs = build_tfim_observables(args.n, args.J, args.h, periodic)
    est_obs = {k: obs[k] for k in ("Mz2", "Mx", "Mzz")}
    exact = exact_tfim_observables(args.n, args.J, args.h, t_final, periodic)
    print(f"[exact] t={t_final:g} {exact}")

    service = get_service()
    backend = None
    binfo: Dict[str, Any] = {}
    if service and not args.no_hardware:
        try:
            backend = service.backend(args.backend) if args.backend \
                else select_backend(service)
        except Exception as ex:
            print(f"[backend] hardware unavailable ({ex}); Aer only.")
            backend = None
    ref_backend = backend or (get_fake_backend(args.fake) if args.fake else None)
    if ref_backend is not None:
        binfo = backend_properties_dict(ref_backend)

    layout_info: Dict[str, Any] = {}
    layout = None
    if ref_backend is not None:
        try:
            _, layout_info = hardware_aware_transpile(core, ref_backend, args.n, periodic)
            layout = layout_info["layout"]
        except Exception as ex:
            print(f"[layout] skipped ({ex})")
    rows = transpile_compare(core, ref_backend, initial_layout=layout)

    noise = noise_model_like(ref_backend) if ref_backend is not None else noise_model_like()
    aer_vals = {str(r["optimization_level"]): observables_from_evs(
        run_expectation_aer(core, est_obs, shots=args.shots, seed=seed, noise_model=noise,
                            optimization_level=r["optimization_level"]))
        for r in rows}
    ideal = observables_from_evs(run_expectation_aer(core, est_obs, shots=args.shots, seed=seed))
    zne = run_with_zne_aer(core, est_obs, shots=args.shots, seed=seed, noise_model=noise)
    counts = run_counts_aer(qc, shots=args.shots, seed=seed, noise_model=noise)
    from_counts = observables_from_counts(counts, args.n, periodic)

    hw: Dict[str, Any] = {}
    if backend and not args.no_hardware:
        try:
            hw["no_mit"] = run_expectation_hardware(core, est_obs, backend, shots=args.shots,
                                                    resilience_level=0, initial_layout=layout)
            hw["zne"] = run_expectation_hardware(core, est_obs, backend, shots=args.shots,
                                                 use_zne=True, initial_layout=layout)
        except Exception as ex:
            print(f"[hw] execution failed: {ex}")
            hw["error"] = str(ex)

    analysis = analyze(rows, exact, aer_vals, zne)
    payload = {
        "meta": {"time": datetime.now(timezone.utc).isoformat(), "seed": seed,
                 "params": vars(args), "t": t_final,
                 "versions": _versions(), "wall_s": round(time.time() - t0, 2)},
        "backend": binfo,
        "transpilation": rows,
        "layout": layout_info,
        "exact": exact,
        "aer_ideal": ideal,
        "aer_noisy_by_level": aer_vals,
        "zne_aer": zne,
        "aer_counts_observables": from_counts,
        "hardware": hw,
        "analysis": analysis,
    }
    save_json(payload, outdir / "tfim_results.json")
    make_plots(analysis, outdir / "figs")
    print(f"[done] {time.time()-t0:.1f}s → {outdir}/")
    return 0


def _versions() -> Dict[str, str]:
    out = {}
    for mod in ("qiskit", "qiskit_aer", "qiskit_ibm_runtime", "numpy", "scipy"):
        try:
            out[mod] = __import__(mod).__version__  # type: ignore
        except Exception:
            try:
                import importlib.metadata as md
                out[mod] = md.version(mod.replace("_", "-"))
            except Exception:
                out[mod] = "unknown"
    return out


if __name__ == "__main__":
    raise SystemExit(main())
