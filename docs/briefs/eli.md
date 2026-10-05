# Eli: Iceberg error-detection simulation

**You own:** `src/qfest/iceberg.py`.
**Read first:** the Iceberg paper (`arXiv:2211.06703`), pages 1–3 and Figure 1 only.

## The idea in plain words
- We store k "logical" qubits inside k + 2 physical qubits. The two extras are called t and b.
- Two parity checks should always give +1: "the product of X on every qubit" (S_X) and "the product of Z on every qubit" (S_Z). If a check gives −1, an error happened, so we **throw that shot away**.
- Two extra helper qubits (ancillas) measure those checks in the middle of the circuit.
- Every gate we need becomes a 2-qubit physical gate:
  - logical ZZ on i, j → `rzz` on physical i, j
  - logical X on i → `rxx` on (t, i)
  - logical Z on i → `rzz` on (b, i)
- **Our question:** how many shots do we throw away, and is the answer better than without the code?

## Steps
1. **`initialisation(k)` (day 1).** Hadamard on t, then a chain of CNOTs to make (|00…0⟩ + |11…1⟩)/√2.
   *Test:* compare against that state with `Statevector`.
2. **`logical_rotation(...)` (day 1–2).** Use the table above, plus the paper's supplementary Eqs. 1–12.
   *Test:* the gate keeps the S_X and S_Z checks at +1.
3. **`syndrome_round(k)` (day 2).** Start simple: one ancilla collects the Z parity and one collects the X parity. Upgrade to the paper's ABBB…BA order only if time allows.
4. **`encode_tfim(...)` and `decode_counts(...)` (day 2).** Build the encoded Trotter circuit, then post-process: drop shots where an ancilla reads 1 or the Z parity is odd.
   *Test:* with no noise, results match `tfim_circuit` and the discard rate is 0.
5. **Noisy sweep (day 3).** Use Toto's noise model. Vary the number of Trotter steps and syndrome rounds (0, 1, 2). Record the discard rate and the Mzz error vs ED, and save JSON.
6. **Hand-off.** Give your circuits to Boutaina for the heavy-hex SWAP analysis.

David pairs with you on step 2. Ask early if stuck for more than 30 minutes.
