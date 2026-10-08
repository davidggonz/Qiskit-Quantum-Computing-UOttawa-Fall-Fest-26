# Results: zero-noise extrapolation of TFIM dynamics on IBM hardware

*Draft for the final write-up (David). Numbers come from `results/zne_time_sweep_ibm_*.json`;
figures from `experiments/plot_zne_summary.py`. References in [brackets] are listed at the end.*

## Summary

On IBM's `ibm_quebec` (Heron, heavy-hex), we ran a 12-spin transverse-field Ising quench
on a natively embedded ring with up to 2400 two-qubit gates. Exponential zero-noise
extrapolation (ZNE) cut the mean error of Mz, Mx and Mzz from 0.11–0.14 (raw) to 0.03, a
3.8–5× reduction. The result held in three independent runs at two field strengths (h/J = 1
and 2). Linear ZNE, the method of the previous project [1], removed only 17–26% of the error
(0.10–0.12) and stayed 7–45σ away from exact diagonalization at t ≥ 1.5. Each sweep used
under a minute of QPU time.

![Summary](../results/figs/zne_summary.png)

## 1. Setup

| Item | Value |
|---|---|
| Model | 1D TFIM, H = −J Σ ZᵢZᵢ₊₁ − h Σ Xᵢ, periodic ring, N = 12, J = 1, h/J ∈ {1, 2} [3] |
| Protocol | Quench from \|0…0⟩, times t = 0.5, 1.0, 1.5, 2.0 |
| Circuit | 2nd-order Trotter [4], Δt = 0.1, edge-coloured RZZ layers, fused Rx half-steps (`qfest.tfim`) |
| Observables | Mz = √⟨(ΣZ)²⟩/N (RMS, as in [1]), Mx = ⟨ΣX⟩/N, Mzz = ⟨ΣZZ⟩/N |
| Reference | Exact diagonalization (`qfest.ed`) |
| Device | `ibm_quebec`, qubits 47–51, 58, 71–67, 57 (a 12-qubit cycle of the heavy-hex lattice [10]) |
| Error suppression | TREX readout mitigation [8], Pauli twirling [6], XpXm dynamical decoupling [7] (Runtime resilience level 1) |
| ZNE | Local unitary folding of every CZ, G → G(G†G)ᵏ [11]; λ = 1–5, even λ by partial folding; 4000 shots per circuit |
| Extrapolators | Linear, Richardson (quadratic) [5], exponential a·e^(−bλ)+c [12], with c free or fixed at the fully depolarized value (Mz → √(1/N), Mx = Mzz → 0) |

**Hardware mapping.** The heavy-hex lattice has no cycle shorter than 12 qubits, so N = 12 is
the smallest periodic chain that runs without SWAP gates. For three Trotter steps the ring
needs the minimum 72 CZ gates. The original layout routine needed 105 (+46%).

## 2. Main result: exponential ZNE is 3.8–5× more accurate than raw, linear is not

Mean |estimate − exact| over Mz, Mx, Mzz and the four times (12 points per run):

| Estimator | h/J = 1, run 1 (λ = 1,3,5) | h/J = 1, run 2 (λ = 1–5) | h/J = 2 (λ = 1–5) |
|---|---|---|---|
| Raw (λ = 1) | 0.142 | 0.145 | 0.114 |
| Linear ZNE [1] | 0.105 | 0.118 | 0.095 |
| Richardson (quadratic) | 0.059 | 0.061 | 0.050 |
| **Exponential, free asymptote** | **0.033** | **0.029** | **0.030** |
| **Exponential, fixed asymptote** | **0.029** | **0.029** | **0.029** |
| Reduction, raw → exponential (free) | 4.4× | 5.0× | 3.8× |

- **Reproducible.** Runs 1 and 2 were independent jobs a day apart, and both give 0.03. Their
  raw λ = 1 values differ by 0.002–0.016, up to 4.1σ of combined shot noise. We therefore quote
  ±0.01 of device drift as a systematic uncertainty on top of the statistical error bars.
- **Linear ZNE is precise but biased.** Its error bars are ±0.003–0.011, yet at t ≥ 1.5 it
  stays 7–45σ from exact (≥ 19σ for Mz and Mzz), and it never lands above exact
  (0 of 36 points). That's the same failure mode the previous report described as
  "begins to overcorrect after t ≈ 1.5" [1], now measured with uncertainties on IBM hardware.
- **Exponential ZNE holds past t = 1.5.** At h/J = 1 (run 2), t = 2.0, it gives
  Mz 0.517 ± 0.027 against exact 0.500, Mx 0.471 ± 0.034 against 0.494, and Mzz 0.471 ± 0.030
  against 0.506. That's in the regime where the λ = 5 circuit has 2400 two-qubit gates.

![Error vs time](../results/figs/zne_summary_vs_time.png)

## 3. Is the exponential model justified? Fit diagnostics

With five noise factors, both exponential fits have spare degrees of freedom, so we can test
them (χ²/dof ≈ 1 means consistent):

| χ²/dof, free asymptote (2 dof) | t = 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|
| h/J = 1: Mz / Mx / Mzz | 16.9 / 0.4 / 3.7 | 0.4 / 0.1 / 1.2 | 2.7 / 1.0 / 0.3 | 0.5 / 0.7 / 0.4 |
| h/J = 2: Mz / Mx / Mzz | 0.0 / 0.3 / 0.4 | 4.4 / 1.0 / 5.7 | 7.2 / 0.0 / 0.8 | 0.2 / 0.4 / 0.2 |

- **Most fits pass.** 19 of 24 have χ²/dof ≤ 2.7, and Mx passes everywhere.
- **The failures are localized:** h/J = 1 at t = 0.5 (Mz, Mzz) and h/J = 2 at t = 1.0–1.5.
  - At the shallowest point (120 two-qubit gates), noise that CZ folding does not amplify
    (single-qubit gates, idling, residual readout error) is a large share of the total. This
    leaves a 0.02–0.03 floor that no extrapolator removes.
  - At h/J = 2, t = 1.0, the raw Mz² drops steeply and then flattens
    (0.228 → 0.172 → 0.144 → 0.115 → 0.108). That pattern suggests two decay rates rather
    than one.
- **Why the free asymptote is the main estimator.** The fixed asymptote assumes fully
  depolarizing noise. For Mz this is slightly wrong: T1 relaxation pulls the noisy limit
  above √(1/N), and the measured Mz² at λ = 5, t = 2 is 0.089, above 1/N = 0.083.

## 4. Where the remaining error comes from

- **Trotter error is small.** At h/J = 2, the noiseless 2nd-order Trotter circuit differs from
  exact diagonalization by at most 0.011 over these times. Measured against the noiseless
  Trotter circuit instead of ED, the exponential-ZNE error is 0.030–0.032. The residual is
  hardware noise, not discretization.
- **Systematic underestimate at h/J = 1.** Raw and linear estimates lie below exact at all
  36 points. At h/J = 1 the exponential estimate lies below exact at 23 of 24 points, the
  signature of noise that folding does not amplify, consistent with the t = 0.5 floor.
  At h/J = 2 the exponential estimate scatters on both sides (4 of 12 above).
- **Mzz at late times at h/J = 2** has the largest residual (0.06–0.10). There Mzz falls to
  0.022 (λ = 4) and 0.006 (λ = 5), almost fully scrambled, so the top noise factors carry
  little information. A smaller λ range (1–3) would suit this regime.

## 5. Cost

| Job | Circuits | λ | QPU time | Resource units |
|---|---|---|---|---|
| `db2qm2m8v0ts73c3c1eg` (h/J = 1, run 1) | 12 | 1, 3, 5 | 39 s | 0.65 |
| `db2r7us2ljfc73d51kig` (h/J = 1, run 2) | 20 | 1–5 | 58 s | 0.96 |
| `db2rctc7f06c73aqdi10` (h/J = 2) | 20 | 1–5 | not yet queried | – |

That works out to about 10 s of QPU time per time point with 3 noise factors and about 15 s
with 5. Run 2 cost 1.49× run 1 rather than the 1.67× its gate counts suggest, because part
of each job is fixed overhead.

## 6. Comparison with the previous project [1]

| Topic | Previous (Quantinuum H2 emulators: H2-1LE, H2-Emulator) | This work (IBM `ibm_quebec`) |
|---|---|---|
| ZNE extrapolation | Linear, λ = 1, 3, 5; "begins to overcorrect after t ≈ 1.5" | Exponential stays within ~0.03 through t = 2, at 2 field strengths, 3 runs |
| Linear ZNE (same circuits) | Main method | Removes 17–26% of error; 7–45σ biased at t ≥ 1.5 |
| Uncertainties | Not reported for TFIM ZNE | Shot-noise error bars, χ² fit tests, run-to-run drift |
| Hardware mapping | All-to-all ions, no routing | Native 12-qubit heavy-hex cycle, 0 SWAPs |
| Table 1 (ED baseline) | Contains transcription errors | Corrected; see `docs/TABLE1_DISCREPANCIES.md` |

The previous report's quantitative ZNE figures (42.3% → 16.7%) are for Fermi–Hubbard, not the
TFIM, so the TFIM comparison is qualitative: where their linear ZNE degraded, ours is replaced
by an extrapolator that does not.

## 7. Limitations

1. One device (`ibm_quebec`), one system size (N = 12), four time points per field.
2. Folding amplifies CZ noise only. Single-qubit, idle and residual readout errors are not
   extrapolated away, which explains the ~0.02–0.03 floor.
3. The exponential model fails its fit test at 5 of 24 points (Section 3), so it's a good
   description of the noise but not an exact one.
4. Device drift between runs (±0.01) is comparable to the late-time error bars.
5. No quantum advantage is claimed: N = 12 is exactly solvable classically, which is what lets
   us check the results.

## 8. Reproducing

```bash
python experiments/zne_time_sweep_ibm.py plan --backend ibm_quebec
python experiments/zne_time_sweep_ibm.py submit --backend ibm_quebec [--h 2.0] [--lambdas 1,3,5]
python experiments/zne_time_sweep_ibm.py analyze --job-id <id> [--subset 1,3,5]
python experiments/zne_time_sweep_ibm.py analyze --from-json results/<file>.json   # no QPU needed
python experiments/plot_zne_summary.py
```

## References

*Check volume and page numbers before they go on slides.*

1. Quantum in Silico, "Simulation of Materials for Next-Generation Energy Devices," Quantathon
   CR 2026, Challenge 3, technical report (2026).
2. C. N. Self, M. Benedetti, D. Amaro, "Protecting expressive circuits with a quantum error
   detection code," Nature Physics (2024), arXiv:2211.06703.
3. P. Pfeuty, "The one-dimensional Ising model with a transverse field," Ann. Phys. 57, 79 (1970).
4. M. Suzuki, "Generalized Trotter's formula and systematic approximants of exponential
   operators," Commun. Math. Phys. 51, 183 (1976).
5. K. Temme, S. Bravyi, J. M. Gambetta, "Error mitigation for short-depth quantum circuits,"
   Phys. Rev. Lett. 119, 180509 (2017).
6. J. J. Wallman, J. Emerson, "Noise tailoring for scalable quantum computation via randomized
   compiling," Phys. Rev. A 94, 052325 (2016).
7. L. Viola, S. Lloyd, "Dynamical suppression of decoherence in two-state quantum systems,"
   Phys. Rev. A 58, 2733 (1998).
8. E. van den Berg, Z. K. Minev, K. Temme, "Model-free readout-error mitigation for quantum
   expectation values," Phys. Rev. A 105, 032620 (2022).
9. Y. Kim et al., "Evidence for the utility of quantum computing before fault tolerance,"
   Nature 618, 500 (2023).
10. C. Chamberland et al., "Topological and subsystem codes on low-degree graphs with flag
    qubits," Phys. Rev. X 10, 011022 (2020).
11. T. Giurgica-Tiron et al., "Digital zero noise extrapolation for quantum error mitigation,"
    IEEE QCE 2020, arXiv:2005.10921.
12. S. Endo, S. C. Benjamin, Y. Li, "Practical quantum error mitigation for near-future
    applications," Phys. Rev. X 8, 031027 (2018).
13. A. Javadi-Abhari et al., "Quantum computing with Qiskit," arXiv:2405.08810 (2024).
