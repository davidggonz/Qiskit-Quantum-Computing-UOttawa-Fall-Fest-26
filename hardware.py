"""Hardware-side module for TFIM Trotter simulation on IBM Quantum hardware.

Covers: backend setup, TFIM circuit generation, transpilation comparison
(0-3), hardware-aware layout via best qubits, ZNE error mitigation,
Aer + hardware execution, exact-diagonalization analysis, JSON + PNG outputs.

Qiskit 1.0+ / 2.x compatible. No hardcoded API keys (.env + env vars).
Reproducible via --seed.

Usage (Colab / local):
    pip install "qiskit[visualization]==2.5.2" "qiskit-aer==0.17.2" \
      "qiskit-ibm-runtime==0.50.0" scipy scikit-learn matplotlib pandas
    cp .env.example .env   # then edit QISKIT_IBM_TOKEN
    python hardware.py --n 4 --steps 3 --shots 4096 --no-hardware
    python hardware.py --n 4 --steps 3 --shots 4096 --backend ibm_marrakesh
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import random
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from qiskit import QuantumCircuit, transpile
    from qiskit.quantum_info import SparsePauliOp
    _QISKIT_AVAILABLE = True
except Exception:  # pragma: no cover - import-time guard for docs/CI
    QuantumCircuit = None  # type: ignore
    SparsePauliOp = None  # type: ignore
    transpile = None  # type: ignore
    _QISKIT_AVAILABLE = False


# ---------------------------------------------------------------------------
# Env / service
# ---------------------------------------------------------------------------

def load_env_file(path: str | Path = ".env") -> None:
    """Load KEY=VALUE lines from a .env file into os.environ (no override)."""
    p = Path(path)
    if not p.exists():
        return
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip().strip('"').strip("'")
        if k and k not in os.environ:
            os.environ[k] = v


def get_api_token(env_vars: Sequence[str] = ("QISKIT_IBM_TOKEN", "PINQ2_TOKEN", "IBM_QUANTUM_TOKEN")) -> str:
    """Return first non-empty token from env vars, else ''."""
    for name in env_vars:
        tok = os.environ.get(name, "").strip().strip('"').strip("'")
        if tok and "PASTE" not in tok and "HERE" not in tok:
            return tok
    return ""


def get_service(channel: str = "ibm_quantum_platform", token: str = "", instance: str | None = None):
    """Connect to IBM Quantum via QiskitRuntimeService (lazy import)."""
    from qiskit_ibm_runtime import QiskitRuntimeService

    load_env_file()
    token = token or get_api_token()
    if not token:
        raise RuntimeError("No API token found. Set QISKIT_IBM_TOKEN in .env or env.")
    kwargs: dict[str, Any] = {"channel": channel, "token": token}
    inst = instance or os.environ.get("QISKIT_IBM_INSTANCE", "").strip()
    if inst:
        kwargs["instance"] = inst
    return QiskitRuntimeService(**kwargs)


def select_backend(service, preferred: Sequence[str] = ("ibm_marrakesh", "ibm_fez"), min_qubits: int = 4):
    """Select first available backend from preferred list, else least-busy.

    Falls back to any operational backend with >= min_qubits.
    """
    backends = service.backends()
    by_name = {getattr(b, "name", str(b)): b for b in backends}
    for name in preferred:
        if name in by_name:
            return by_name[name]
    cands = [b for b in backends if getattr(b, "num_qubits", 0) >= min_qubits]
    if not cands:
        return backends[0]
    try:
        return min(cands, key=lambda b: service.jobs(limit=5, backend_name=getattr(b, "name", "")).count() if hasattr(service, "jobs") else 0)
    except Exception:
        return cands[0]


def backend_info_dict(backend) -> dict[str, Any]:
    """Collect backend properties into a JSON-serializable dict."""
    info: dict[str, Any] = {}
    try:
        info["name"] = getattr(backend, "name", str(backend))
    except Exception:
        info["name"] = str(backend)
    try:
        info["num_qubits"] = int(backend.num_qubits)
    except Exception:
        pass
    try:
        cm = backend.coupling_map
        info["coupling_map"] = [list(e) for e in cm.get_edges()] if cm is not None else None
    except Exception:
        info["coupling_map"] = None
    try:
        info["basis_gates"] = sorted(list(backend.target.operation_names))  # type: ignore
    except Exception:
        try:
            info["basis_gates"] = list(backend.configuration().basis_gates)
        except Exception:
            info["basis_gates"] = None
    try:
        info["dt"] = getattr(backend, "dt", None)
    except Exception:
        pass
    # Error rates (best-effort; properties() needs network).
    try:
        props = backend.properties()
        rq, g1, g2 = [], {}, {}
        for i, q in enumerate(props.qubits):
            try:
                ro = next(p for p in q if p.name == "readout_error").value
                rq.append(float(ro))
            except Exception:
                pass
        info["readout_error_mean"] = float(np.mean(rq)) if rq else None
        for g in props.gates:
            try:
                err = next(p for p in g.parameters if p.name == "gate_error").value
                key = (g.gate, tuple(g.qubits))
                if len(g.qubits) == 1:
                    g1[key] = float(err)
                elif len(g.qubits) == 2:
                    g2[key] = float(err)
            except Exception:
                continue
        if g1:
            info["1q_error_mean"] = float(np.mean(list(g1.values())))
        if g2:
            info["2q_error_mean"] = float(np.mean(list(g2.values())))
            info["2q_error_min"] = float(min(g2.values()))
            info["2q_error_max"] = float(max(g2.values()))
    except Exception as ex:
        info["properties_error"] = f"{type(ex).__name__}"
    return info


def print_backend_properties(backend) -> dict[str, Any]:
    """Print backend properties and return them as dict."""
    info = backend_info_dict(backend)
    print(f"Backend: {info.get('name')}")
    print(f"  num_qubits: {info.get('num_qubits')}")
    print(f"  basis_gates: {info.get('basis_gates')}")
    print(f"  coupling_map edges: {len(info.get('coupling_map') or [])}")
    for k in ("readout_error_mean", "1q_error_mean", "2q_error_mean", "2q_error_min", "2q_error_max"):
        if info.get(k) is not None:
            print(f"  {k}: {info[k]:.5f}")
    return info


# ---------------------------------------------------------------------------
# TFIM circuit + observables
# ---------------------------------------------------------------------------

def build_tfim_circuit(n: int = 4, J: float = 1.0, h: float = 0.5,
                       dt: float = 0.1, steps: int = 3,
                       add_measure: bool = True, barrier: bool = True,
                       init_plus: bool = False) -> Any:
    """Build first-order Trotter circuit for 1D TFIM.

    H = -J * sum_<i,j> Z_i Z_j - h * sum_i X_i, open chain.
    U(dt) ~= prod exp(+i J dt ZZ) prod exp(+i h dt X).
    RZZ(phi)=exp(-i phi/2 ZZ) -> phi=-2*J*dt; RX(phi)=exp(-i phi/2 X) -> phi=-2*h*dt.
    Init: |0>^N by default (ferromagnetic reference, matches demo baseline).
    Set init_plus=True for H on all qubits (|+>^N, large-h reference).
    """
    if not _QISKIT_AVAILABLE:
        raise RuntimeError("qiskit is required")
    if n not in (4, 5) and not (2 <= n <= 12):
        raise ValueError("N should be 4 or 5 for the workshop (2..12 allowed).")
    qc = QuantumCircuit(n, n if add_measure else 0)
    if init_plus:
        qc.h(range(n))
        if barrier:
            qc.barrier()
    theta_zz = -2.0 * J * dt
    theta_x = -2.0 * h * dt
    for _ in range(steps):
        for i in range(n - 1):
            qc.rzz(theta_zz, i, i + 1)
        for i in range(n):
            qc.rx(theta_x, i)
        if barrier:
            qc.barrier()
    if add_measure:
        qc.measure(range(n), range(n))
    return qc


def build_tfim_circuit_2nd_order(n, J, h, dt, steps, periodic=True, initial_state="quench"):
    """Build a symmetric second-order TFIM circuit with fused X half-rotations."""
    if not _QISKIT_AVAILABLE:
        raise RuntimeError("qiskit is required")
    if n % 2:
        raise ValueError("N must be even for canonical TFIM + Iceberg compatibility")
    if n < 2:
        raise ValueError("n must be at least 2")
    if initial_state not in ("quench", "adiabatic"):
        raise ValueError("initial_state must be 'quench' or 'adiabatic'")

    qc = QuantumCircuit(n)
    if initial_state == "adiabatic":
        qc.h(range(n))

    edges = [(i, i + 1) for i in range(n - 1)]
    if periodic and n > 2:
        edges.append((n - 1, 0))
    edge_layers = [[], []]
    for index, edge in enumerate(edges):
        edge_layers[index % 2].append(edge)

    theta_x_half = -h * dt
    theta_x_full = -2.0 * h * dt
    theta_zz = -2.0 * J * dt
    qc.rx(theta_x_half, range(n))
    for _ in range(steps):
        for layer in edge_layers:
            for left, right in layer:
                qc.rzz(theta_zz, left, right)
        qc.rx(theta_x_half if _ == steps - 1 else theta_x_full, range(n))
    return qc


def build_tfim_circuit_1st_order(n, J, h, dt, steps, periodic=False):
    """Legacy first-order builder, kept for the appendix ablation."""
    if not _QISKIT_AVAILABLE:
        raise RuntimeError("qiskit is required")
    if n < 2:
        raise ValueError("n must be at least 2")
    qc = QuantumCircuit(n)
    edges = [(i, i + 1) for i in range(n - 1)]
    if periodic and n > 2:
        edges.append((n - 1, 0))
    for _ in range(steps):
        for left, right in edges:
            qc.rzz(-2.0 * J * dt, left, right)
        qc.rx(-2.0 * h * dt, range(n))
    return qc


def measure_observables(qc, n):
    """Return measurement circuits for M_z, M_x, and nearest-neighbor M_zz."""
    if not _QISKIT_AVAILABLE:
        raise RuntimeError("qiskit is required")
    if qc.num_qubits != n:
        raise ValueError("n must match the circuit's qubit count")
    z_circuit = qc.copy()
    z_circuit.measure_all()
    x_circuit = qc.copy()
    x_circuit.h(range(n))
    x_circuit.measure_all()
    return {"M_z": z_circuit, "M_x": x_circuit, "M_zz": z_circuit}


def tfim_observables(n: int) -> dict[str, Any]:
    """Return magnetization + ZZ correlator observables (SparsePauliOp)."""
    if not _QISKIT_AVAILABLE:
        raise RuntimeError("qiskit is required")
    z_terms = []
    for i in range(n):
        s = ["I"] * n
        s[n - 1 - i] = "Z"  # Qiskit little-endian: qubit i -> string pos n-1-i
        z_terms.append(("".join(s), 1.0 / n))
    zz_terms = []
    for i in range(n - 1):
        s = ["I"] * n
        s[n - 1 - i] = "Z"
        s[n - 1 - (i + 1)] = "Z"
        zz_terms.append(("".join(s), 1.0 / max(1, n - 1)))
    return {
        "magnetization": SparsePauliOp.from_list(z_terms),
        "zz_correlator": SparsePauliOp.from_list(zz_terms),
    }


# ---------------------------------------------------------------------------
# Transpilation
# ---------------------------------------------------------------------------

CX_LIKE = {"cx", "cnot", "cz", "ecr", "rzz", "crx", "cry", "crz", "zx", "iswap"}


def _stats_of(circ) -> dict[str, Any]:
    ops = dict(circ.count_ops())
    return {
        "depth": int(circ.depth()),
        "qubits": int(circ.num_qubits),
        "total_gates": int(circ.size()),
        "cx_like": int(sum(v for k, v in ops.items() if k.lower() in CX_LIKE)),
        "counts": {str(k): int(v) for k, v in ops.items()},
    }


def estimate_circuit_error(circ, backend=None) -> float | None:
    """Rough 1 - prod(1-err) estimate using backend properties. None if unavailable."""
    if backend is None:
        return None
    try:
        props = backend.properties()
    except Exception:
        return None
    try:
        errmap: dict[tuple, float] = {}
        for g in props.gates:
            try:
                e = next(p for p in g.parameters if p.name == "gate_error").value
                errmap[(g.gate, tuple(g.qubits))] = float(e)
            except Exception:
                continue
        layout = getattr(circ, "_layout", None)
        # Without exact mapping, use mean per gate type.
        from collections import defaultdict
        by_type: dict[str, list[float]] = defaultdict(list)
        for (gate, _), e in errmap.items():
            by_type[gate].append(e)
        mean_err = {g: float(np.mean(v)) for g, v in by_type.items()}
        fail = 1.0
        for gate, count in circ.count_ops().items():
            e = mean_err.get(gate, mean_err.get(gate.lower(), 0.0))
            fail *= (1.0 - e) ** int(count)
        return float(1.0 - fail)
    except Exception:
        return None


def transpile_comparison(circuit, backend=None, levels: Sequence[int] = (0, 1, 2, 3),
                         seed: int = 42, vary_seed: bool = False) -> list[dict[str, Any]]:
    """Transpile at each optimization level; record depth/gates/CNOT/est error."""
    rows = []
    for lvl in levels:
        kw: dict[str, Any] = {
            "optimization_level": int(lvl),
            "seed_transpiler": seed + int(lvl) if vary_seed else seed,
        }
        if backend is not None:
            kw["backend"] = backend
        t = transpile(circuit, **kw)
        row = {"optimization_level": int(lvl)}
        row.update(_stats_of(t))
        row["est_error"] = estimate_circuit_error(t, backend)
        rows.append(row)
    return rows


def print_comparison_table(rows: Sequence[dict[str, Any]]) -> None:
    print(f"{'lvl':>3} | {'depth':>5} | {'gates':>5} | {'cx-like':>7} | {'est_err':>9}")
    print("-" * 44)
    for r in rows:
        e = r.get("est_error")
        print(f"{r['optimization_level']:>3} | {r['depth']:>5} | {r['total_gates']:>5} | "
              f"{r['cx_like']:>7} | {(f'{e:.4f}' if e is not None else 'n/a'):>9}")


def find_best_qubits(backend, n: int) -> list[int]:
    """Pick N connected physical qubits with lowest readout+1q+2q error.

    Greedy: rank qubits by readout+1q error, BFS-expand over coupling graph
    picking the lowest-cost connected set. Falls back to range(n).
    """
    try:
        props = backend.properties()
        cm = backend.coupling_map
        edges = [tuple(e) for e in cm.get_edges()]
    except Exception:
        return list(range(n))
    try:
        adj: dict[int, set[int]] = {}
        for a, b in edges:
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
        qcost: dict[int, float] = {}
        for i in range(backend.num_qubits):
            try:
                q = props.qubit(i)
                ro = next(p for p in q if p.name == "readout_error").value
            except Exception:
                ro = 0.02
            qcost[i] = float(ro)
        ecost: dict[tuple[int, int], float] = {}
        for g in props.gates:
            if len(g.qubits) == 2:
                try:
                    e = next(p for p in g.parameters if p.name == "gate_error").value
                    a, b = int(g.qubits[0]), int(g.qubits[1])
                    ecost[(a, b)] = ecost.get((a, b), float(e))
                except Exception:
                    continue

        def set_cost(nodes: set[int]) -> float:
            c = sum(qcost.get(q, 0.02) for q in nodes)
            # add internal edge costs (approx: MST weight via min edges)
            internal = [ecost.get((a, b), ecost.get((b, a), 0.02))
                        for a in nodes for b in nodes if a < b and (b in adj.get(a, set()))]
            internal.sort()
            c += sum(internal[: max(0, len(nodes) - 1)])
            return c

        best: set[int] | None = None
        best_c = float("inf")
        for start in sorted(qcost, key=qcost.get)[: min(12, len(qcost))]:
            nodes = {start}
            while len(nodes) < n:
                frontier = {w for q in nodes for w in adj.get(q, set())} - nodes
                if not frontier:
                    break
                nxt = min(frontier, key=lambda w: set_cost(nodes | {w}))
                nodes.add(nxt)
            if len(nodes) == n and set_cost(nodes) < best_c:
                best, best_c = set(nodes), set_cost(nodes)
        if best:
            # order along a chain for linear TFIM
            ordered = [next(iter(best))]
            rem = set(best) - set(ordered)
            while rem:
                last = ordered[-1]
                nxts = [w for w in adj.get(last, set()) if w in rem]
                nxt = min(nxts, key=lambda w: ecost.get((last, w), ecost.get((w, last), 0.02))) if nxts else min(rem)
                ordered.append(nxt)
                rem.remove(nxt)
            return [int(q) for q in ordered]
    except Exception:
        pass
    return list(range(n))


def hardware_aware_transpile(circuit, backend, n: int, seed: int = 42) -> tuple[Any, dict[str, Any]]:
    """Transpile opt-level 3 with initial_layout=best qubits; return (circ, meta)."""
    best = find_best_qubits(backend, n)
    t = transpile(circuit, backend=backend, optimization_level=3,
                  initial_layout=best, seed_transpiler=seed)
    meta = {"best_qubits": [int(q) for q in best]}
    meta.update(_stats_of(t))
    meta["est_error"] = estimate_circuit_error(t, backend)
    return t, meta


# ---------------------------------------------------------------------------
# ZNE (manual folding; works on Aer + hardware via estimator values)
# ---------------------------------------------------------------------------

def fold_circuit_global(circuit, scale: int):
    """Global unitary folding: U -> U (U^dagger U)^((scale-1)/2). Odd scales only."""
    if scale == 1:
        return circuit.copy()
    if scale % 2 == 0:
        raise ValueError("scale must be odd (1, 3, 5, ...)")
    base = circuit.copy()
    try:
        inv = circuit.inverse()
    except Exception:
        inv = circuit.copy().inverse()
    out = base.copy()
    for _ in range((scale - 1) // 2):
        out = out.compose(inv).compose(base)
    return out


def zne_extrapolate(noise_factors: Sequence[float], values: Sequence[float],
                    method: str = "linear") -> float:
    """Extrapolate to zero noise. linear (polyfit deg1) or exponential."""
    x = np.asarray(noise_factors, dtype=float)
    y = np.asarray(values, dtype=float)
    if method == "linear" or len(x) < 3:
        p = np.polyfit(x, y, 1)
        return float(p[1])
    # exponential: y = a + b*exp(-c*x) via scipy
    try:
        from scipy.optimize import curve_fit

        def f(xx, a, b, c):
            return a + b * np.exp(-c * xx)

        popt, _ = curve_fit(f, x, y, p0=(y[-1], y[0] - y[-1], 1.0), maxfev=5000)
        return float(popt[0])
    except Exception:
        p = np.polyfit(x, y, 1)
        return float(p[1])


def _extract_estimator_value(res) -> float:
    for attr in ("data",):
        pass
    try:  # Aer / runtime EstimatorV2: result[0].data.evs
        return float(res.data.evs[0])
    except Exception:
        pass
    try:
        return float(res.values[0])
    except Exception:
        pass
    try:
        return float(np.asarray(res).ravel()[0])
    except Exception as ex:
        raise RuntimeError(f"Cannot extract estimator value: {ex}")


def aer_estimator_value(circuit_no_measure, observable, shots: int = 4096, seed: int = 42,
                        noise_model=None) -> float:
    """Expectation via Aer EstimatorV2 (fallback to statevector)."""
    try:
        from qiskit_aer.primitives import EstimatorV2
        opts: dict[str, Any] = {"shots": shots, "seed": seed}
        est = EstimatorV2(options=opts)
        if noise_model is not None:
            try:
                est = EstimatorV2(options={**opts, "noise_model": noise_model})
            except Exception:
                pass
        job = est.run([(circuit_no_measure, observable)])
        return _extract_estimator_value(job.result()[0])
    except Exception:
        from qiskit.quantum_info import Statevector
        sv = Statevector(circuit_no_measure.remove_final_measurements(inplace=False))
        return float(sv.expectation_value(observable).real)


def run_zne_aer(circuit_no_measure, observable, shots: int = 4096, seed: int = 42,
                noise_factors: Sequence[int] = (1, 3, 5), extrapolator: str = "linear",
                noise_model=None) -> dict[str, Any]:
    """Manual ZNE on Aer: fold circuit at each scale, estimate, extrapolate."""
    raw = []
    for s in noise_factors:
        folded = fold_circuit_global(circuit_no_measure, int(s))
        raw.append(aer_estimator_value(folded, observable, shots=shots, seed=seed, noise_model=noise_model))
    return {"noise_factors": list(noise_factors), "raw": [float(v) for v in raw],
            "mitigated": float(zne_extrapolate(noise_factors, raw, extrapolator)),
            "extrapolator": extrapolator}


# ---------------------------------------------------------------------------
# Execution: Aer counts + hardware
# ---------------------------------------------------------------------------

def counts_magnetization(counts: dict[str, int], n: int) -> float:
    """<Z_avg> from counts (bitstring -> +1 for '0', -1 for '1')."""
    shots = sum(counts.values())
    tot = 0.0
    for bitstr, c in counts.items():
        s = bitstr.replace(" ", "")
        z_sum = sum(1.0 if b == "0" else -1.0 for b in s[-n:])
        tot += (z_sum / n) * c
    return float(tot / shots) if shots else 0.0


def counts_observables(counts_z: dict[str, int], counts_x: dict[str, int],
                       n: int, periodic: bool = True) -> dict[str, float]:
    """Compute M_z, M_x, and nearest-neighbor M_zz from measurement counts."""
    shots = sum(counts_z.values())
    if not shots:
        return {"M_z": 0.0, "M_x": 0.0, "M_zz": 0.0}
    edge_count = n if periodic else max(1, n - 1)
    zz_total = 0.0
    for bitstring, count in counts_z.items():
        bits = bitstring.replace(" ", "")[-n:]
        z_values = [1 if bits[-1 - i] == "0" else -1 for i in range(n)]
        zz_total += count * sum(
            z_values[i] * z_values[(i + 1) % n]
            for i in range(edge_count)
        ) / edge_count
    return {
        "M_z": counts_magnetization(counts_z, n),
        "M_x": counts_magnetization(counts_x, n),
        "M_zz": float(zz_total / shots),
    }


def run_aer_counts(circuit_with_measure, shots: int = 4096, seed: int = 42) -> dict[str, int]:
    from qiskit_aer import AerSimulator
    sim = AerSimulator(seed_simulator=seed)
    res = sim.run(circuit_with_measure, shots=shots, seed_simulator=seed).result()
    return dict(res.get_counts())


def _get_fake_backend(name: str):
    from qiskit_ibm_runtime.fake_provider import FakeFez, FakeMarrakesh

    return {"ibm_marrakesh": FakeMarrakesh, "ibm_fez": FakeFez}[name]()


def _get_aer_simulator(fake_backend, noise_mode: str):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel
    from src.qfest.noise import from_backend, simple_model

    if noise_mode == "backend":
        if fake_backend is None:
            raise ValueError("--noise-model backend requires --fake-backend")
        noise_model = from_backend(fake_backend)
    elif noise_mode == "simple":
        noise_model = simple_model()
    else:
        noise_model = NoiseModel()
    if fake_backend is not None:
        return AerSimulator.from_backend(fake_backend, noise_model=noise_model)
    return AerSimulator(noise_model=noise_model)


def _run_zne_counts(circuit, simulator, backend, n: int, shots: int, seed: int,
                    periodic: bool, noise_factors: Sequence[int] = (1, 3, 5)) -> dict[str, Any]:
    """Run unchanged global unitary folds through a sampled Aer execution path."""
    raw = []
    for index, scale in enumerate(noise_factors):
        folded = fold_circuit_global(circuit, int(scale))
        measured = folded.copy()
        measured.measure_all()
        transpiled = transpile(
            measured, backend=backend, optimization_level=3,
            seed_transpiler=seed + index,
        )
        counts = dict(simulator.run(
            transpiled, shots=shots, seed_simulator=seed + index
        ).result().get_counts())
        raw.append(counts_observables(counts, counts, n, periodic)["M_z"])
    return {
        "noise_factors": list(noise_factors),
        "raw": [float(value) for value in raw],
        "mitigated": float(zne_extrapolate(noise_factors, raw, "linear")),
        "extrapolator": "linear",
    }


def run_hardware(service, backend_name: str, circuit_no_measure, observables: dict[str, Any],
                 shots: int = 4096, use_zne_option: bool = False) -> dict[str, Any]:
    """Run on real backend via EstimatorV2 (+ Sampler counts). Returns expvals."""
    from qiskit_ibm_runtime import EstimatorV2, SamplerV2, Session
    backend = service.backend(backend_name)
    out: dict[str, Any] = {"backend": getattr(backend, "name", backend_name)}
    with Session(backend=backend) as session:
        # Estimator
        try:
            from qiskit_ibm_runtime.options import EstimatorOptions
            opts = EstimatorOptions()
            opts.default_shots = shots
            if use_zne_option:
                try:
                    opts.resilience_level = 2
                    opts.resilience.zne_mitigation = True
                except Exception:
                    pass
            est = EstimatorV2(mode=session, options=opts)
        except Exception:
            est = EstimatorV2(mode=session)
        pubs = [(circuit_no_measure, [obs]) for obs in observables.values()]
        # EstimatorV2.run accepts list of pubs; flatten per observable for clarity
        expvals: dict[str, float] = {}
        for key, obs in observables.items():
            job = est.run([(circuit_no_measure, obs)])
            expvals[key] = _extract_estimator_value(job.result()[0])
        out["expectations"] = expvals
        # Counts via Sampler (needs measurements)
        try:
            meas = circuit_no_measure.copy()
            meas.measure_all()
            sampler = SamplerV2(mode=session)
            job = sampler.run([meas], shots=shots)
            res = job.result()[0]
            try:
                counts = res.data.meas.get_counts()
            except Exception:
                counts = res.data["meas"].get_counts()
            if isinstance(counts, list):
                counts = counts[0]
            out["counts"] = {str(k): int(v) for k, v in dict(counts).items()}
        except Exception as ex:
            out["counts_error"] = f"{type(ex).__name__}: {ex}"
    return out


# ---------------------------------------------------------------------------
# Exact baseline
# ---------------------------------------------------------------------------

def exact_magnetization(n: int, J: float, h: float, dt: float, steps: int,
                        init_plus: bool = False, periodic: bool = False) -> float:
    """Exact <Z_avg> via numpy expm, same init as circuit (default |0>^N)."""
    from scipy.linalg import expm
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I = np.eye(2, dtype=complex)

    def kron_n(ops):
        m = ops[0]
        for o in ops[1:]:
            m = np.kron(m, o)
        return m

    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(n - 1):
        ops = [I] * n
        ops[i], ops[i + 1] = Z, Z
        H += -J * kron_n(ops)
    if periodic and n > 2:
        ops = [I] * n
        ops[n - 1], ops[0] = Z, Z
        H += -J * kron_n(ops)
    for i in range(n):
        ops = [I] * n
        ops[i] = X
        H += -h * kron_n(ops)
    U = expm(-1j * H * dt * steps)
    if init_plus:
        plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
        psi0 = plus
        for _ in range(n - 1):
            psi0 = np.kron(psi0, plus)
    else:
        zero = np.array([1, 0], dtype=complex)
        psi0 = zero
        for _ in range(n - 1):
            psi0 = np.kron(psi0, zero)
    psi = U @ psi0
    mag = 0.0
    for i in range(n):
        ops = [I] * n
        ops[i] = Z
        Zop = kron_n(ops)
        mag += float(np.vdot(psi, Zop @ psi).real)
    return float(mag / n)


# ---------------------------------------------------------------------------
# Analysis / plots / save
# ---------------------------------------------------------------------------

def ensure_outdir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_json(payload: dict[str, Any], path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, default=str))
    return p


def make_plots(analysis_rows: list[dict], zne: dict | None, exact: float, outdir: str | Path) -> list[str]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    lvls = [r["optimization_level"] for r in analysis_rows if isinstance(r.get("optimization_level"), int)]
    acc = [abs(float(r.get("aer_value", 0)) - exact) for r in analysis_rows if isinstance(r.get("optimization_level"), int)]
    depths = [r["depth"] for r in analysis_rows if isinstance(r.get("optimization_level"), int)]

    fig, ax = plt.subplots()
    ax.plot(lvls, acc, marker="o")
    ax.set_xlabel("optimization_level")
    ax.set_ylabel("|Aer - exact| magnetization")
    ax.set_title("Accuracy vs optimization level")
    ax.grid(True, alpha=0.3)
    p1 = out / "accuracy_vs_level.png"
    fig.savefig(p1, dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.bar([str(l) for l in lvls], depths)
    ax.set_xlabel("optimization_level")
    ax.set_ylabel("depth")
    ax.set_title("Depth vs optimization level")
    p2 = out / "depth_vs_level.png"
    fig.savefig(p2, dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()
    if zne:
        x = list(zne["noise_factors"])
        y = list(zne["raw"])
        ax.plot(x, y, marker="o", label="raw")
        ax.axhline(zne["mitigated"], linestyle="--", label="zero-noise extrapolation")
        ax.scatter([0], [zne["mitigated"]], marker="*", s=140, zorder=4,
                   label="extrapolated at zero noise")
        spread = max(y) - min(y)
        padding = max(spread * 0.35, 0.002)
        ax.set_ylim(min(min(y), zne["mitigated"]) - padding,
                    max(max(y), zne["mitigated"]) + padding)
        ax.set_xlabel("noise factor")
        ax.set_ylabel("M_z estimate")
        ax.set_title(f"ZNE M_z (exact={exact:.4f})")
        ax.legend()
        ax.grid(True, alpha=0.3)
    p3 = out / "error_vs_mitigation.png"
    fig.savefig(p3, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [str(p1), str(p2), str(p3)]


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(n: int = 4, J: float = 1.0, h: float = 0.5, dt: float = 0.1,
                 steps: int = 3, shots: int = 4096, seed: int = 42,
                 backend_name: str = "", channel: str = "ibm_quantum_platform",
                 outdir: str = "results", no_hardware: bool = False,
                 with_noise_demo: bool = False, init_plus: bool = False,
                 trotter_order: int = 2, periodic: bool = True,
                 protocol: str = "quench", fake_backend_name: str = "",
                 noise_mode: str = "ideal") -> dict[str, Any]:
    """Full pipeline; hardware section skipped if no token / --no-hardware."""
    t0 = time.time()
    random.seed(seed)
    np.random.seed(seed)
    out = ensure_outdir(outdir)

    if trotter_order == 2:
        circ = build_tfim_circuit_2nd_order(n, J, h, dt, steps, periodic, protocol)
    else:
        circ = build_tfim_circuit_1st_order(n, J, h, dt, steps, periodic)
        if protocol == "adiabatic":
            initial = QuantumCircuit(n)
            initial.h(range(n))
            circ = initial.compose(circ)
    circ_meas = circ.copy()
    circ_meas.measure_all()
    obs = tfim_observables(n)

    backend, binfo = None, {}
    fake_backend = _get_fake_backend(fake_backend_name) if fake_backend_name else None
    if not no_hardware and fake_backend is None:
        try:
            svc = get_service(channel=channel)
            backend = select_backend(svc, min_qubits=n) if not backend_name else svc.backend(backend_name)
            binfo = print_backend_properties(backend)
        except Exception as ex:
            print(f"Hardware setup skipped: {type(ex).__name__}: {ex}")
            backend = None
    # Transpilation (needs backend or basis-agnostic)
    target_backend = fake_backend or backend
    rows = transpile_comparison(
        circ_meas, target_backend, levels=(0, 1, 2, 3), seed=seed,
        vary_seed=fake_backend is not None,
    )
    print_comparison_table(rows)

    layout: dict[str, Any] = {}
    if backend is not None:
        try:
            _, layout = hardware_aware_transpile(circ_meas, backend, n, seed=seed)
            print(f"Best qubits: {layout.get('best_qubits')} depth={layout.get('depth')}")
        except Exception as ex:
            layout = {"error": f"{type(ex).__name__}: {ex}"}

    exact = exact_magnetization(n, J, h, dt, steps,
                                init_plus=(protocol == "adiabatic"), periodic=periodic)
    noisy_path = fake_backend is not None or noise_mode != "ideal"
    simulator = _get_aer_simulator(fake_backend or backend, noise_mode) if noisy_path else None
    measured = measure_observables(circ, n)
    noisy_observables: dict[str, float] = {}
    noisy_counts: dict[str, int] = {}
    if noisy_path:
        aer_vals = {}
        for row in rows:
            level = int(row["optimization_level"])
            level_seed = seed + level
            z_circuit = transpile(
                measured["M_z"], backend=target_backend,
                optimization_level=level, seed_transpiler=level_seed,
            )
            x_circuit = transpile(
                measured["M_x"], backend=target_backend,
                optimization_level=level, seed_transpiler=level_seed,
            )
            z_counts = dict(simulator.run(
                z_circuit, shots=shots, seed_simulator=level_seed
            ).result().get_counts())
            x_counts = dict(simulator.run(
                x_circuit, shots=shots, seed_simulator=level_seed
            ).result().get_counts())
            level_observables = counts_observables(z_counts, x_counts, n, periodic)
            aer_vals[str(level)] = level_observables["M_z"]
            if level == 3:
                noisy_observables = level_observables
                noisy_counts = z_counts
    else:
        aer_vals = {str(r["optimization_level"]): aer_estimator_value(circ, obs["magnetization"], shots, seed) for r in rows}

    noise_model = None
    if with_noise_demo:
        try:
            from qiskit_aer.noise import NoiseModel, depolarizing_error
            nm = NoiseModel()
            nm.add_all_qubit_quantum_error(depolarizing_error(0.01, 1), ["rx", "x", "h"])
            nm.add_all_qubit_quantum_error(depolarizing_error(0.03, 2), ["cx", "rzz", "ecr"])
            noise_model = nm
        except Exception:
            noise_model = None
    if noisy_path:
        zne = _run_zne_counts(
            circ, simulator, target_backend, n, shots, seed, periodic, (1, 3, 5)
        )
        counts = noisy_counts
    else:
        zne = run_zne_aer(circ, obs["magnetization"], shots, seed, (1, 3, 5), "linear", noise_model)
        counts = run_aer_counts(circ_meas, shots, seed)
    mag_counts = counts_magnetization(counts, n)
    if noisy_observables:
        print("Noisy observables (optimization level 3): " +
              ", ".join(f"{key}={value:.6f}" for key, value in noisy_observables.items()))
    if fake_backend is not None:
        ideal_value = aer_estimator_value(circ, obs["magnetization"], shots, seed)
        noisy_value = aer_vals["3"]
        delta = abs(noisy_value - ideal_value)
        print(f"Fake-backend noise check: ideal={ideal_value:.6f} noisy={noisy_value:.6f} delta={delta:.6g}")
        if noise_mode != "ideal" and delta <= 1e-6:
            print("WARNING: noisy raw expectation matches ideal within 1e-6; stopping.")
            raise RuntimeError("Fake-backend noise was not applied to the execution path")

    hw: dict[str, Any] = {}
    if backend is not None:
        try:
            import qiskit_ibm_runtime as _rt  # noqa: F401
            svc = get_service(channel=channel)
            bname = getattr(backend, "name", backend_name)
            hw = run_hardware(svc, bname, circ, {"magnetization": obs["magnetization"]}, shots, use_zne_option=False)
        except Exception as ex:
            hw = {"error": f"{type(ex).__name__}: {ex}"}

    table = []
    for r in rows:
        v = aer_vals[str(r["optimization_level"])]
        table.append({**r, "aer_value": float(v), "abs_error": float(abs(v - exact))})
    table.append({"optimization_level": "zne", "mitigated": zne["mitigated"], "raw": zne["raw"]})

    try:
        import qiskit as _q, qiskit_aer as _a, qiskit_ibm_runtime as _r
        versions = {"qiskit": _q.__version__, "qiskit_aer": _a.__version__,
                    "qiskit_ibm_runtime": _r.__version__, "numpy": np.__version__}
        try:
            import scipy as _s
            versions["scipy"] = _s.__version__
        except Exception:
            pass
    except Exception:
        versions = {}
    payload = {
        "meta": {"time": _dt.datetime.now(_dt.timezone.utc).isoformat(), "seed": seed,
                 "params": {"n": n, "J": J, "h": h, "dt": dt, "steps": steps, "shots": shots,
                            "trotter_order": trotter_order, "periodic": periodic, "protocol": protocol,
                            "backend": fake_backend_name or (getattr(backend, "name", backend_name) if backend is not None else backend_name),
                            "no_hardware": backend is None, "outdir": str(outdir),
                            "init_plus": init_plus},
                 "versions": versions, "wall_s": round(time.time() - t0, 2)},
        "backend": binfo, "transpilation": rows, "layout": layout,
        "exact_magnetization": exact, "aer_expectations": aer_vals,
        "zne_aer": zne, "aer_counts_mag": mag_counts, "hardware": hw,
        "analysis": {"exact": exact, "table": table},
    }
    save_json(payload, out / "tfim_results.json")
    plots = make_plots(table, zne, exact, out)
    print(f"Saved JSON + {len(plots)} plots to {out}")
    print(f"Exact mag={exact:.4f} Aer={list(aer_vals.values())[0]:.4f} ZNE={zne['mitigated']:.4f} counts-mag={mag_counts:.4f}")
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="TFIM hardware module")
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--J", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=None, help="Directly set h instead of using --h-over-j")
    ap.add_argument("--h-over-j", type=float, choices=(0.5, 1.0, 2.0), default=0.5)
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--steps", type=int, default=20)
    ap.add_argument("--trotter-order", type=int, choices=(1, 2), default=2)
    topology = ap.add_mutually_exclusive_group()
    topology.add_argument("--periodic", action="store_true", dest="periodic")
    topology.add_argument("--open-chain", action="store_false", dest="periodic")
    ap.set_defaults(periodic=True)
    ap.add_argument("--protocol", choices=("quench", "adiabatic"), default="quench")
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--backend", type=str, default="")
    ap.add_argument("--channel", type=str, default="ibm_quantum_platform")
    ap.add_argument("--outdir", type=str, default="results")
    ap.add_argument("--no-hardware", action="store_true")
    ap.add_argument("--with-noise-demo", action="store_true")
    ap.add_argument("--fake-backend", choices=("none", "ibm_marrakesh", "ibm_fez"), default="none")
    ap.add_argument("--noise-model", "--noise", dest="noise_model",
                    choices=("ideal", "simple", "backend"), default=None)
    ap.add_argument("--init-plus", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args(argv)
    if args.n % 2:
        raise ValueError("N must be even for canonical TFIM + Iceberg compatibility")
    h = args.h if args.h is not None else args.J * args.h_over_j
    protocol = "adiabatic" if args.init_plus else args.protocol
    fake_backend_name = "" if args.fake_backend == "none" else args.fake_backend
    noise_mode = args.noise_model or ("backend" if fake_backend_name else "ideal")
    run_pipeline(args.n, args.J, h, args.dt, args.steps, args.shots, args.seed,
                 args.backend, args.channel, args.outdir, args.no_hardware, args.with_noise_demo,
                 protocol == "adiabatic", args.trotter_order, args.periodic, protocol,
                 fake_backend_name, noise_mode)
    return 0


if __name__ == "__main__":
    sys.exit(main())
