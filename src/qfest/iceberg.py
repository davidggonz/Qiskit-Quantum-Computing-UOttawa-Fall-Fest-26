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


def initialisation(k):
    raise NotImplementedError


def syndrome_round(k):
    raise NotImplementedError


def logical_rotation(kind, i, j, theta):
    raise NotImplementedError


def encode_tfim(n, J, h, dt, steps, rounds):
    raise NotImplementedError


def decode_counts(counts, k):
    raise NotImplementedError
