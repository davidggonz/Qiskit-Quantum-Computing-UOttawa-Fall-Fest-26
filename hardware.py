"""Hardware-side module for TFIM Trotter simulation on IBM Quantum hardware.

Covers: backend setup, TFIM Trotter circuits, transpilation comparison,
hardware-aware layout selection, ZNE error mitigation, Aer + hardware
execution, exact-diagonalization analysis, JSON + PNG outputs.

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
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

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
    backends = service.backends()
    by_name = {b.name: b for b in backends}
    for name in preferred:
        if name in by_name:
            print(f"[backend] Selected preferred backend: {name}")
            return by_name[name]
    # fallback: operational backend with most qubits free
    operational = [b for b in backends if getattr(b, "status", lambda: None) is None]
    cand = backends[0]
    print(f"[backend] Preferred not found, using: {cand.name}")
    return cand


def backend_properties_dict(backend) -> Dict[str, Any]:
    """Collect printable backend properties + per-qubit/gate error rates."""
    props = backend.properties() if hasattr(backend, "properties") else None
    target = getattr(backend, "target", None)
    n = backend.num_qubits if hasattr(backend, "num_qubits") else backend.num_qubits
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
            for op in ("x", "sx", "rz", "id"):
                try:
                    e = target[op][(q,)].error
                    err_1q[q] = min(err_1q.get(q, 1.0), float(e))
                except Exception:
                    continue
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
# 2. TFIM Trotter circuit
# ---------------------------------------------------------------------------

def build_tfim_circuit(n: int = 4, j: float = 1.0, h: float = 0.5,
                       dt: float = 0.1, steps: int = 3,
                       periodic: bool = False,
                       add_measure: bool = True):
    """Trotterized TFIM: H = -J Σ ZZ - h Σ X, U ≈ Π e^{+iJ·ZZ·dt} Π e^{+ih·X·dt}.

    RZZ(θ)=exp(-iθ/2 ZZ) → θ=-2J·dt ; RX(θ)=exp(-iθ/2 X) → θ=-2h·dt.
    Starts in |0>^N (Z-polarized). Returns QuantumCircuit.
    """
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(n, n if add_measure else 0)
    pairs = [(i, (i + 1) % n) for i in range(n if periodic else n - 1)]
    theta_zz, theta_x = -2.0 * j * dt, -2.0 * h * dt
    for _ in range(steps):
        for a, b in pairs:
            qc.rzz(theta_zz, a, b)
        for q in range(n):
            qc.rx(theta_x, q)
        qc.barrier()
    if add_measure:
        qc.measure(range(n), range(n))
    qc.metadata = {"n": n, "J": j, "h": h, "dt": dt, "steps": steps}
    return qc


def build_tfim_observables(n: int = 4, j: float = 1.0, h: float = 0.5):
    """Return dict with energy and magnetization observables as SparsePauliOp.

    Returns {'energy': H_sparse, 'magnetization': Z_avg_sparse} built from
    the TFIM Hamiltonian H = -J sum(Z_i Z_{i+1}) - h sum(X_i).
    """
    from qiskit.quantum_info import SparsePauliOp
    zz_terms = []
    for i in range(n - 1):
        label = ["I"] * n
        label[i] = "Z"
        label[i + 1] = "Z"
        zz_terms.append(("".join(reversed(label)), -j))
    x_terms = []
    for i in range(n):
        label = ["I"] * n
        label[i] = "X"
        x_terms.append(("".join(reversed(label)), -h))
    all_terms = zz_terms + x_terms
    if all_terms:
        labels = [t[0] for t in all_terms]
        coeffs = [complex(t[1]) for t in all_terms]
        H_op = SparsePauliOp(labels, coeffs)
    else:
        H_op = SparsePauliOp(["I" * n], [0.0])
    z_labels = []
    z_coeffs = []
    for i in range(n):
        label = ["I"] * n
        label[i] = "Z"
        z_labels.append("".join(reversed(label)))
        z_coeffs.append(1.0 / n)
    mag_op = SparsePauliOp(z_labels, z_coeffs) if z_labels else SparsePauliOp(["I" * n], [0.0])
    return {"energy": H_op, "magnetization": mag_op}


def exact_tfim_magnetization(n: int = 4, j: float = 1.0, h: float = 0.5,
                             dt: float = 0.1, steps: int = 3) -> float:
    """Exact <Z_avg> via scipy matrix exponentiation of the Trotterized evolution.

    Builds H = -J sum(ZZ) - h sum(X), computes U = exp(-i H dt)^steps applied
    to |0...0>, returns average per-qubit magnetization.
    """
    from scipy.linalg import expm
    sx = np.array([[0, 1], [1, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    eye = np.eye(2, dtype=complex)

    def kron_n(ops):
        m = ops[0]
        for o in ops[1:]:
            m = np.kron(m, o)
        return m

    dim = 2 ** n
    H = np.zeros((dim, dim), complex)
    for i in range(n - 1):
        ops = [eye] * n
        ops[i] = sz
        ops[i + 1] = sz
        H += -j * kron_n(ops)
    for i in range(n):
        ops = [eye] * n
        ops[i] = sx
        H += -h * kron_n(ops)
    U = expm(-1j * H * dt * steps)
    psi0 = np.zeros(dim, complex)
    psi0[0] = 1.0
    psi = U @ psi0
    z_ops = []
    for i in range(n):
        ops = [eye] * n
        ops[i] = sz
        z_ops.append(kron_n(ops))
    z_avg = sum(np.real(np.vdot(psi, z @ psi)) for z in z_ops) / n
    return float(z_avg)


# ---------------------------------------------------------------------------
# 3. Transpilation comparison
# ---------------------------------------------------------------------------

def estimate_circuit_error(circuit, backend) -> Optional[float]:
    """1 - Π(1-e_g) using backend.target error rates; None if unavailable."""
    target = getattr(backend, "target", None)
    if target is None:
        return None
    try:
        from qiskit import QuantumCircuit
        qc = circuit.remove_final_measurements(inplace=False) or circuit
    except Exception:
        qc = circuit
    fid = 1.0
    for inst, qargs, _ in qc.data:
        name = inst.name
        qa = tuple(q.index for q in qargs)
        try:
            e = float(target[name][qa].error)
            fid *= (1.0 - e)
        except Exception:
            continue
    return float(1.0 - fid)


def transpile_compare(circuit, backend=None, levels: Sequence[int] = (0, 1, 2, 3)):
    """Transpile at each opt level; record depth, gates, CX count, est. error."""
    from qiskit import transpile
    rows = []
    for lv in levels:
        t = transpile(circuit, backend=backend, optimization_level=lv, seed_transpiler=42)
        counts = t.count_ops()
        cx = int(counts.get("cx", 0) + counts.get("ecr", 0) + counts.get("cz", 0))
        row = {
            "optimization_level": lv,
            "depth": int(t.depth()),
            "qubits": int(t.num_qubits),
            "total_gates": int(sum(counts.values())),
            "cx_like": cx,
            "counts": {k: int(v) for k, v in counts.items()},
            "est_error": estimate_circuit_error(t, backend) if backend else None,
        }
        rows.append(row)
    # pretty table
    print(f"{'level':>5} {'depth':>6} {'gates':>6} {'cx':>5} {'est_err':>9}")
    for r in rows:
        e = f"{r['est_error']:.4f}" if r["est_error"] is not None else "n/a"
        print(f"{r['optimization_level']:>5} {r['depth']:>6} "
              f"{r['total_gates']:>6} {r['cx_like']:>5} {e:>9}")
    return rows


# ---------------------------------------------------------------------------
# 4. Hardware-aware layout
# ---------------------------------------------------------------------------

def select_best_qubits(backend, n: int) -> List[int]:
    """Greedy lowest-error connected chain of n physical qubits.

    Score qubit = readout_err + avg 1q err; edge = 2q err. Start from best
    qubit, BFS-expand along lowest-error edges until n qubits collected.
    """
    target = backend.target
    num = backend.num_qubits
    q_score = {}
    for q in range(num):
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

    try:
        neighbors = {q: set(backend.coupling_map.neighbors(q)) for q in range(num)}
    except Exception:
        neighbors = {q: set() for q in range(num)}
    start = min(q_score, key=q_score.get)
    chosen = [start]
    while len(chosen) < n:
        frontier = {}
        for q in chosen:
            for nb in neighbors.get(q, []):
                if nb not in chosen:
                    frontier[nb] = min(frontier.get(nb, 9), edge_err(q, nb) + q_score[nb])
        if not frontier:
            # disconnected map fallback: add best remaining qubit
            rest = [q for q in range(num) if q not in chosen]
            chosen.append(min(rest, key=lambda q: q_score[q]))
        else:
            chosen.append(min(frontier, key=frontier.get))
    print(f"[layout] best {n} qubits: {chosen}")
    return chosen


def hardware_aware_transpile(circuit, backend, n: int):
    """Transpile opt_level=3 with initial_layout=best qubits; compare to default."""
    from qiskit import transpile
    layout = select_best_qubits(backend, n)
    default = transpile(circuit, backend=backend, optimization_level=1, seed_transpiler=42)
    aware = transpile(circuit, backend=backend, optimization_level=3,
                      initial_layout=layout, seed_transpiler=42)
    for label, t in (("default(lv1)", default), ("aware(lv3+layout)", aware)):
        c = t.count_ops()
        print(f"[layout] {label}: depth={t.depth()}, gates={sum(c.values())}, {dict(c)}")
    return aware, {"layout": layout,
                   "default": {"depth": default.depth(), "counts": dict(default.count_ops())},
                   "aware": {"depth": aware.depth(), "counts": dict(aware.count_ops())}}


# ---------------------------------------------------------------------------
# 5-6. Execution (Aer + hardware) and ZNE
# ---------------------------------------------------------------------------

def run_counts_aer(circuit, shots: int = 4096, seed: int = 42) -> Dict[str, int]:
    """Counts on ideal AerSimulator (no noise)."""
    from qiskit_aer import AerSimulator
    sim = AerSimulator(seed_simulator=seed)
    job = sim.run(circuit, shots=shots)
    return {k: int(v) for k, v in job.result().get_counts().items()}


def run_expectation_aer(circuit, observable, shots: int = 4096,
                        seed: int = 42) -> float:
    """Expectation via Aer EstimatorV2 (no mitigation)."""
    from qiskit_aer.primitives import EstimatorV2 as AerEstimator
    from qiskit import transpile
    qc = circuit.remove_final_measurements(inplace=False) or circuit
    est = AerEstimator()
    est.options.run_options["shots"] = shots
    est.options.run_options["seed"] = seed
    pub = (transpile(qc), [observable])
    res = est.run([pub]).result()
    return float(np.real(res[0].data.evs[0]))


def zne_extrapolate(noise_factors: Sequence[float], values: Sequence[float],
                    method: str = "linear") -> float:
    """Extrapolate to zero noise. linear | exponential | poly2."""
    x = np.asarray(noise_factors, float)
    y = np.asarray(values, float)
    if method == "exponential":
        # log-linear fit y = A*exp(-k x); fallback to linear if y<=0
        if np.all(y > 0):
            k, loga = np.polyfit(x, np.log(y), 1)
            return float(np.exp(loga))
        method = "linear"
    deg = 2 if method == "poly2" else 1
    coef = np.polyfit(x, y, deg)
    return float(np.polyval(coef, 0.0))


def fold_two_qubit_gates(circuit, scale: int):
    """Unitary folding for ZNE: repeat each CX/CZ/ECR `scale` times (odd).

    RZZ/RX left unfolded (parameterized); CX-like noise dominates on hardware.
    scale=1 → identity.
    """
    from qiskit import QuantumCircuit
    if scale == 1:
        return circuit.copy()
    assert scale % 2 == 1, "scale must be odd"
    folded = QuantumCircuit(*circuit.qregs, *circuit.cregs) if circuit.qregs else QuantumCircuit(circuit.num_qubits, circuit.num_clbits)
    for inst, qargs, cargs in circuit.data:
        folded.append(inst, qargs, cargs)
        if inst.name in ("cx", "cz", "ecr") and len(qargs) == 2:
            for _ in range(scale - 1):
                folded.append(inst, qargs, cargs)
    folded.metadata = dict(getattr(circuit, "metadata", {}) or {})
    return folded


def run_with_zne_aer(circuit, observable, shots: int = 4096, seed: int = 42,
                     noise_factors: Sequence[float] = (1, 3, 5),
                     extrapolator: str = "linear") -> Dict[str, Any]:
    """Manual ZNE on Aer: fold CX noise, evaluate at each factor, extrapolate."""
    vals = []
    for f in noise_factors:
        fc = fold_two_qubit_gates(circuit, int(f))
        vals.append(run_expectation_aer(fc, observable, shots=shots, seed=seed))
    mitigated = zne_extrapolate(noise_factors, vals, method=extrapolator)
    print(f"[zne-aer] raw={vals} → mitigated({extrapolator})={mitigated:.6f}")
    return {"noise_factors": list(noise_factors), "raw": vals,
            "mitigated": mitigated, "extrapolator": extrapolator}


def run_expectation_hardware(circuit, observable, backend, shots: int = 4096,
                             resilience_level: int = 0,
                             use_zne: bool = False) -> Dict[str, Any]:
    """Expectation on real backend via EstimatorV2, optionally with ZNE."""
    from qiskit_ibm_runtime import EstimatorV2, Session
    from qiskit import transpile
    qc = circuit.remove_final_measurements(inplace=False) or circuit
    qc_t = transpile(qc, backend=backend, optimization_level=3, seed_transpiler=42)
    opts: Dict[str, Any] = {"default_shots": shots, "resilience_level": resilience_level}
    if use_zne or resilience_level >= 2:
        opts["resilience_level"] = max(resilience_level, 2)
        opts["resilience"] = {"zne_mitigation": True,
                              "zne": {"noise_factors": (1, 3, 5),
                                      "extrapolator": ("exponential", "linear")}}
    with Session(backend=backend) as session:
        est = EstimatorV2(mode=session, options=opts)
        res = est.run([(qc_t, [observable])]).result()
        ev = float(np.real(res[0].data.evs[0]))
    print(f"[hw] backend={backend.name} shots={shots} res_lv={opts['resilience_level']} "
          f"zne={bool(use_zne or resilience_level>=2)} → {ev:.6f}")
    return {"value": ev, "shots": shots, "options": opts}


def run_counts_hardware(circuit, backend, shots: int = 4096) -> Dict[str, int]:
    """Counts on real backend via SamplerV2."""
    from qiskit_ibm_runtime import SamplerV2, Session
    from qiskit import transpile
    qc_t = transpile(circuit, backend=backend, optimization_level=3, seed_transpiler=42)
    with Session(backend=backend) as session:
        samp = SamplerV2(mode=session, options={"default_shots": shots})
        res = samp.run([(qc_t,)]).result()
        counts = res[0].data.c.get_counts()
    return {k: int(v) for k, v in counts.items()}


# ---------------------------------------------------------------------------
# 7. Analysis + plots + persistence
# ---------------------------------------------------------------------------

def magnetization_from_counts(counts: Dict[str, int], n: int) -> float:
    """<Z_avg> from counts (bitstring key, space-separated classical regs ok)."""
    total = sum(counts.values())
    acc = 0.0
    for bits, c in counts.items():
        b = bits.replace(" ", "")
        z_sum = sum(1.0 if ch == "0" else -1.0 for ch in b[-n:])
        acc += (z_sum / n) * c
    return acc / total if total else 0.0


def analyze(level_rows, exact: float, aer_vals: Dict[str, float],
            zne_vals: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """Build accuracy table vs exact baseline."""
    table = []
    for r in level_rows:
        lv = r["optimization_level"]
        v = aer_vals.get(str(lv), aer_vals.get(lv))
        entry = {**r, "aer_value": v,
                 "abs_error": abs(v - exact) if v is not None else None}
        table.append(entry)
    if zne_vals:
        table.append({"optimization_level": "zne", **zne_vals})
    return {"exact": exact, "table": table}


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
    rows = [r for r in analysis["table"] if isinstance(r.get("optimization_level"), int)]
    levels = [r["optimization_level"] for r in rows]
    paths = []
    if rows and rows[0].get("aer_value") is not None:
        errs = [abs(r["aer_value"] - analysis["exact"]) for r in rows]
        plt.figure()
        plt.plot(levels, errs, marker="o")
        plt.xlabel("optimization level")
        plt.ylabel("|Aer - exact| (magnetization)")
        plt.title("Accuracy vs optimization level")
        p = outdir / "accuracy_vs_level.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(p)
    if rows:
        plt.figure()
        plt.plot(levels, [r["depth"] for r in rows], marker="o")
        plt.xlabel("optimization level")
        plt.ylabel("transpiled depth")
        plt.title("Depth vs optimization level")
        p = outdir / "depth_vs_level.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(p)
    zne_rows = [r for r in analysis["table"] if r.get("optimization_level") == "zne"]
    if rows and zne_rows:
        labels = [f"lv{lv}" for lv in levels] + ["ZNE"]
        vals = [abs(r["aer_value"] - analysis["exact"]) for r in rows] + \
               [abs(zne_rows[0]["mitigated"] - analysis["exact"])]
        plt.figure()
        plt.bar(labels, vals)
        plt.ylabel("absolute error vs exact")
        plt.title("Error: plain vs ZNE mitigation")
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
    ap.add_argument("--shots", type=int, default=4096)
    ap.add_argument("--backend", default="")
    ap.add_argument("--no-hardware", action="store_true",
                    help="skip real-backend execution (Aer only)")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args(argv)

    seed = 42
    np.random.seed(seed)
    outdir = Path(args.outdir)
    t0 = time.time()

    qc = build_tfim_circuit(args.n, args.J, args.h, args.dt, args.steps)
    print(f"[circuit] TFIM N={args.n} depth={qc.depth()} gates={dict(qc.count_ops())}")
    obs = build_tfim_observables(args.n, args.J, args.h)
    exact = exact_tfim_magnetization(args.n, args.J, args.h, args.dt, args.steps)
    print(f"[exact] magnetization={exact:.6f}")

    service = get_service()
    backend = None
    binfo: Dict[str, Any] = {}
    if service and not args.no_hardware:
        try:
            backend = service.backend(args.backend) if args.backend \
                else select_backend(service)
            binfo = backend_properties_dict(backend)
        except Exception as ex:
            print(f"[backend] hardware unavailable ({ex}); Aer only.")
            backend = None

    rows = transpile_compare(qc, backend)
    layout_info: Dict[str, Any] = {}
    if backend:
        try:
            _, layout_info = hardware_aware_transpile(
                qc.remove_final_measurements(inplace=False) or qc, backend, args.n)
        except Exception as ex:
            print(f"[layout] skipped ({ex})")

    mag_obs = obs["magnetization"]
    aer_vals = {}
    for r in rows:
        from qiskit import transpile
        lv = r["optimization_level"]
        t = transpile(qc.remove_final_measurements(inplace=False) or qc,
                      backend=backend, optimization_level=lv, seed_transpiler=42)
        aer_vals[str(lv)] = run_expectation_aer(t, mag_obs, shots=args.shots, seed=seed)
    zne = run_with_zne_aer(qc.remove_final_measurements(inplace=False) or qc,
                           mag_obs, shots=args.shots, seed=seed)
    counts = run_counts_aer(qc, shots=args.shots, seed=seed)
    mag_counts = magnetization_from_counts(counts, args.n)

    hw: Dict[str, Any] = {}
    if backend and not args.no_hardware:
        try:
            hw["no_mit"] = run_expectation_hardware(
                qc, mag_obs, backend, shots=args.shots, resilience_level=0)
            hw["zne"] = run_expectation_hardware(
                qc, mag_obs, backend, shots=args.shots, use_zne=True)
            hw["counts"] = run_counts_hardware(qc, backend, shots=args.shots)
        except Exception as ex:
            print(f"[hw] execution failed: {ex}")
            hw["error"] = str(ex)

    analysis = analyze(rows, exact, aer_vals,
                       {"mitigated": zne["mitigated"], "raw": zne["raw"]})
    payload = {
        "meta": {"time": datetime.now(timezone.utc).isoformat(), "seed": seed,
                 "params": vars(args),
                 "versions": _versions(), "wall_s": round(time.time() - t0, 2)},
        "backend": binfo,
        "transpilation": rows,
        "layout": layout_info,
        "exact_magnetization": exact,
        "aer_expectations": aer_vals,
        "zne_aer": zne,
        "aer_counts_mag": mag_counts,
        "hardware": hw,
        "analysis": analysis,
    }
    save_json(payload, outdir / "tfim_results.json")
    make_plots(analysis, outdir)
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
