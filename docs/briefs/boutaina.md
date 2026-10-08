# Boutaina: hardware lead

**You own:** `src/qfest/hardware.py`. That covers backend setup, layout, transpilation, mitigation settings and job submission.
**Time budget: about 7 hours total.** Steps are front-loaded. If you run short, hand the remaining step to David.

## Steps
1. **Account and backend (day 1, 1 h).** Log in to IBM Quantum, check our plan's QPU minutes, and pick a Heron backend. Save your token locally with `QiskitRuntimeService.save_account`. **Never commit the token.**
2. **`get_backend()` (day 1, 1 h).** Return the real backend, or a fake Heron backend from `qiskit_ibm_runtime.fake_provider` when offline, so the team can test without QPU time.
3. **`embed_ring()` and `transpile_report()` (day 1, 1.5 h).** Find a cycle of N physical qubits in the coupling map (rustworkx or networkx cycle search) and use it as `initial_layout`. Transpile `tfim_circuit(...)` and report 2-qubit gate count, depth and SWAP count.
   *Done when:* a 12-qubit ring transpiles with **0 SWAPs**.
4. **Mitigation options (day 2, 1.5 h).** Set up EstimatorV2 with Pauli twirling, dynamical decoupling, TREX readout mitigation, and ZNE with linear and exponential extrapolators. Test on the fake backend. Check whether native RZZ (fractional) gates work with these options; if not, we use CZ/ECR.
5. **`run_batch()` and the hardware run (day 3, 1 h, with David).** Submit everything as **one** batch. Save the job IDs and the raw results to `results/`.
6. **Iceberg vs heavy-hex (day 3–4, 1 h).** When Eli's circuits are ready, transpile them to all-to-all and to heavy-hex and record the extra SWAPs and depth. Give the numbers to Toto for the plot.
