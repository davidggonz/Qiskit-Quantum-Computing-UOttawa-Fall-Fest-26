# Team plan: hardware-aware TFIM simulation on IBM Quantum

Goal, approach and risks: see `proposal/Proposal_EN.pdf` (EN) and `proposal/Propuesta_ES.pdf` (ES).
Time budget: 3-4 days. Core target: 1D TFIM, N = 6-12, Trotter vs multi-product formulas, better
mitigation than linear ZNE, a few IBM hardware runs. Stretch: Iceberg-vs-topology study, Hubbard in simulation.

## Roles

| Person | Role | Module(s) | Tasks |
|---|---|---|---|
| Leader | Physics lead, integrator | `tfim.py`, `mpf.py`, final write-up | Interface contract; Trotter circuits validated vs ED; MPF; review all PRs; analysis and write-up; Hubbard (stretch) |
| Boutaina | Hardware lead (about 7 h/week, front-loaded) | `hardware.py` | Account and backend choice; ring embedding on heavy-hex; transpilation reports; mitigation options (twirling, DD, TREX, ZNE); one batched job submission; Iceberg SWAP/depth analysis |
| Ririsha | Classical baseline and numerics | `ed.py`, `extrapolation.py`, `metrics.py` | ED baseline (done, verify and extend); extrapolators (linear, Richardson, exponential) with tests; error-vs-time metrics |
| Eli | Iceberg simulation | `iceberg.py` | Encoder, syndrome rounds, logical rotations, decoding; Aer noisy sweeps for discard rate and fidelity vs depth |
| Toto | Noise models, plots, docs | `noise.py`, `plots.py`, README | Aer noise model from backend calibration; all figures and tables; reproducibility instructions; slides |

Boutaina's hours are limited: her tasks have a hard deadline at the end of day 1 and the Leader
shadows her as backup. Her day-3 job submission is done together with the Leader.

## Schedule

- **Day 1:** interface contract (done in this skeleton); ED and Trotter circuits validated; Boutaina: account, backend, first transpilation (2-3 h); others: IBM Quantum Learning basics.
- **Day 2:** MPF vs 2nd-order Trotter; extrapolators; Iceberg circuits and sims; noise model; mitigation config review; freeze hardware circuits.
- **Day 3:** batched hardware runs (Boutaina + Leader); simulations continue.
- **Day 4:** plots, comparison table vs the previous report, write-up. No new features.

## Interface contract

- `qfest.tfim.tfim_circuit(n, J, h, dt, steps, order=2, periodic=True) -> QuantumCircuit` (quench from |0...0>, no measurements).
- `qfest.ed.quench(n, times, J, h)` and `ground_state_observables(n, J, h)` return Mz, Mx, Mzz as in the report (Mz is the RMS magnitude).
- Qubit i is the i-th least significant bit (Qiskit order).
- Every experiment saves one JSON via `qfest.results.save`, schema in `src/qfest/results.py`.
- Extrapolators: `f(lams, vals) -> float` in `qfest.extrapolation.EXTRAPOLATORS`.

## Working rules

- One branch per person (`name/topic`), pull requests into `main`, the Leader reviews.
- Every module gets a test in `tests/`; run `pytest` before pushing.
- Never commit IBM API tokens. Use `QiskitRuntimeService.save_account` locally or an environment variable.
- Hardware time is scarce: test everything on Aer or a fake backend first, then batch.
- Known issue: the previous report's Table 1 has typos, see `docs/TABLE1_DISCREPANCIES.md`.

## References to acknowledge

Verify volume and page details before putting them on slides.

- C. N. Self, M. Benedetti, D. Amaro, "Protecting expressive circuits with a quantum error detection code," arXiv:2211.06703 (Nature Physics, 2024).
- Quantum in Silico, "Simulation of Materials for Next-Generation Energy Devices," Quantathon CR 2026, Challenge 3.
- P. Pfeuty, Ann. Phys. 57, 79 (1970). J. Hubbard, Proc. R. Soc. A 276, 238 (1963). P. Jordan, E. Wigner, Z. Phys. 47, 631 (1928).
- M. Suzuki, Commun. Math. Phys. 51, 183 (1976). S. Lloyd, Science 273, 1073 (1996).
- A. Carrera Vazquez et al., "Well-conditioned multi-product formulas for hardware-friendly Hamiltonian simulation," Quantum 7, 1067 (2023).
- Zhuk, Robertson, Bravyi, "Trotter error bounds and dynamic multi-product formulas for Hamiltonian simulation," arXiv:2306.12569.
- Y. Kim et al., "Evidence for the utility of quantum computing before fault tolerance," Nature 618, 500 (2023).
- K. Temme, S. Bravyi, J. Gambetta, PRL 119, 180509 (2017). Y. Li, S. Benjamin, PRX 7, 021050 (2017). Kandala et al., Nature 567, 491 (2019).
- J. Wallman, J. Emerson, PRA 94, 052325 (2016) (twirling). L. Viola, S. Lloyd, PRA 58, 2733 (1998) (dynamical decoupling). van den Berg, Minev, Temme, PRA 105, 032620 (2022) (TREX).
- C. Chamberland et al., PRX 10, 011022 (2020) (heavy-hex, flag qubits).
- A. Javadi-Abhari et al., "Quantum computing with Qiskit," arXiv:2405.08810.
- Stretch: Robledo-Moreno et al., arXiv:2405.05068 (SQD).
