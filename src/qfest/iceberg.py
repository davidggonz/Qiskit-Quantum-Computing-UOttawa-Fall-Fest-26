"""[[k+2, k, 2]] Iceberg code circuits.  (Owner: Eli; Boutaina for transpilation analysis)

Reference: Self, Benedetti, Amaro, arXiv:2211.06703, Fig. 1 and Supplementary Eqs. (1)-(12).
Qubits: data 1..k, plus t, b; ancillas a1 (S_Z flag), a2 (S_X flag).
Logical X_i = X_t X_i, logical Z_i = Z_b Z_i, S_X = X^{k+2}, S_Z = Z^{k+2}.

TODO (in order):
  1. initialisation(k): GHZ preparation of |0bar>^k  (Fig. 1c).
  2. syndrome_round(k): ABBB...BA CNOT pattern with ancillas a1, a2 (Fig. 1d).
  3. logical_rotation(kind, i, j, theta): physical gate from the Eq. (1)-(12) table.
  4. encode_tfim(n, J, h, dt, steps, rounds): encoded Trotter circuit with syndrome rounds.
  5. decode_counts(counts, k): discard on ancilla = 1 or S_Z = -1, then read logical Z_i.
Sanity test: noiseless encoded circuit must give the same observables as tfim.tfim_circuit.
"""

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from .tfim import tfim_circuit

def initialisation(k: int) -> QuantumCircuit:
    if k < 2 or k % 2:
        raise ValueError("k must be even and >= 2")
    t, b = k, k + 1
    qc = QuantumCircuit(k + 2)
    qc.h(t)
    for q in [b, *range(k)]:
        qc.cx(t, q)
    return qc


def syndrome_round(qc: QuantumCircuit, k: int, checks: ClassicalRegister, r: int) -> None:
    a1, a2 = k + 2, k + 3
    for q in range(k + 2):
        qc.cx(q, a1)
    qc.h(a2)
    for q in range(k + 2):
        qc.cx(a2, q)
    qc.h(a2)
    qc.measure(a1, checks[2 * r])
    qc.measure(a2, checks[2 * r + 1])
    qc.reset(a1)
    qc.reset(a2)


def logical_rotation(qc: QuantumCircuit, k: int, kind: str, theta: float, i: int, j: int | None = None) -> None:
    t = k
    if kind == "x":
        qc.rxx(theta, t, i)
    elif kind == "zz":
        qc.rzz(theta, i, j)
    else:
        raise ValueError("kind must be 'x' or 'zz'")


def encode_tfim(n: int, J: float, h: float, dt: float, steps: int, rounds: int = 0) -> QuantumCircuit:
    plain = tfim_circuit(n, J, h, dt, steps)
    qc = initialisation(n)
    if rounds > 0:
        qc.add_register(QuantumRegister(2, "a"))
        checks = ClassicalRegister(2 * rounds, "checks")
        qc.add_register(checks)
    spots = []
    for r in range(rounds):
        spots.append(round((r + 1) * len(plain.data) / (rounds + 1)))
    kinds = {"rx": "x", "rzz": "zz"}
    for idx, ins in enumerate(plain.data):
        qubits = []
        for q in ins.qubits:
            qubits.append(plain.find_bit(q).index)
        logical_rotation(qc, n, kinds[ins.operation.name], ins.operation.params[0], *qubits)
        for r, spot in enumerate(spots):
            if spot == idx + 1:
                syndrome_round(qc, n, checks, r)
    meas = ClassicalRegister(n + 2, "meas")
    qc.add_register(meas)
    qc.measure(range(n + 2), meas)
    return qc


def decode_counts(counts: dict[str, int], k: int) -> tuple[dict[str, int], float]:
    good_runs = {}
    thrown = 0
    for key, val in counts.items():
        meas, *checks = key.split()
        alarm = "1" in "".join(checks)
        if meas.count("1") % 2 == 0 and not alarm:
            b = meas[0]
            magnets = ""
            for c in meas[2:]:
                if b == c:
                    magnets += "0"
                else:
                    magnets += "1"
            good_runs[magnets] = good_runs.get(magnets, 0) + val
        else:
            thrown += val
    return good_runs, thrown / sum(counts.values())
