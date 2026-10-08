# UO-Qubit

**Zero-noise extrapolation of quantum magnetism on a real IBM quantum computer**

## Team

- 🇨🇷 **David Granados** (Universidad de Costa Rica): Team Leader & Physics lead
- 🇲🇦 **Boutaina El Hourri** (ENSA Khouribga): Hardware lead
- 🇱🇧 **Toufic Haddad** (University of Ottawa): Noise models and Plots
- 🇨🇦 **Eli Levasseur** (University of Ottawa): Iceberg Circuit Simulation
- 🇲🇺 **Harshini Gungah**: Classical baseline & Numerics

## Introduction and objective

Designing better energy materials, such as superconductors for efficient power grids,
battery materials and CO₂-capture catalysts, depends on predicting how many interacting
spins and electrons behave. Classical computers hit a wall quickly: in the previous
challenge, exact diagonalization took 0.004 s for 4 spins and 17.5 s for 12, and was no
longer feasible beyond that. Quantum computers are the natural tool for this problem,
but only if their noise can be controlled.

This project builds on the Quantathon CR 2026 Challenge 3 report, *Simulation of
Materials for Next-Generation Energy Devices*. That work used the one-dimensional
transverse-field Ising model (TFIM) as the canonical testbed for magnetic materials:

$$H = -J\sum_i Z_i Z_{i+1} - h\sum_i X_i$$

The first term favours aligned neighbouring spins (ferromagnetic order). The transverse
field adds quantum fluctuations that destroy that order, and the competition produces a
quantum phase transition at h/J = 1. The previous project simulated the time evolution
of this model on Quantinuum's H2 emulators (H2-1LE and H2-Emulator), validated every
circuit against exact diagonalization, and extended the workflow to the 2D Fermi–Hubbard
model, a minimal model of high-temperature superconductors.

The report also documented clear weak spots:

- Its linear zero-noise extrapolation (ZNE, noise factors λ = 1, 3, 5) recovered the exact
  result only until t ≈ 1.5, then began to overcorrect.
- Iceberg error detection discarded up to 95% of the shots (17–42% with minimal checks).
- The mitigated results had no reported error bars.

**Objective:** obtain accurate, statistically validated TFIM dynamics on a real IBM
quantum computer, with an error-mitigation method that keeps working at long evolution
times.

## Our approach

We kept the physics and changed the way it is executed and mitigated:

1. **Real hardware instead of emulators.** All results come from IBM's `ibm_quebec`
   (Heron processor, heavy-hex lattice).
2. **A ring that fits the chip.** The heavy-hex lattice has no closed loop shorter than
   12 qubits, so a 12-spin periodic ring maps onto it with zero SWAP gates. Our layout
   search picks the best of the 21 such rings on the device. For 3 Trotter steps the ring
   needs 72 two-qubit gates, the minimum possible, against 105 with the original greedy
   layout.
3. **A compact circuit.** Second-order Trotter steps (Δt = 0.1) with edge-coloured ZZ
   layers and fused X half-rotations. The circuits run from a quench out of |0…0⟩ to
   t = 0.5–2.0 at h/J = 1 and 2. We measure Mz (RMS, as in the report), Mx and Mzz.
4. **Layered error mitigation.**
   - IBM Runtime error suppression: TREX readout mitigation, Pauli twirling and dynamical
     decoupling.
   - ZNE by folding every CZ gate, CZ → CZ (CZ†CZ)ᵏ, at five noise levels (λ = 1–5).
     Even λ values come from partial folding.
5. **An extrapolation that matches the physics of the noise.** Twirling turns gate errors
   into Pauli noise, which shrinks the signal by a fixed factor per gate. We therefore fit
   ⟨O⟩(λ) ≈ a·e^(−bλ) + c instead of a straight line, and compare it with linear and
   Richardson extrapolation on the same data.
6. **Statistics the previous work lacked.** We report shot-noise error bars, χ² goodness-of-fit
   tests for every extrapolation, and run-to-run reproducibility checks.

## Key results

![Summary of the IBM hardware runs](results/figs/zne_summary.png)

### 1. Exponential ZNE reduces the error 3.8–5× on real hardware

Mean |estimate − exact| over Mz, Mx and Mzz at t = 0.5, 1.0, 1.5 and 2.0 (12 points per run):

| Estimator | h/J = 1, run 1 | h/J = 1, run 2 | h/J = 2 |
|---|---|---|---|
| Raw (no ZNE) | 0.142 | 0.145 | 0.114 |
| Linear ZNE (previous method) | 0.105 | 0.118 | 0.095 |
| **Exponential ZNE (ours)** | **0.033** | **0.029** | **0.030** |
| Reduction, raw → exponential | 4.4× | 5.0× | 3.8× |

**Conclusion:** across three independent runs at two field strengths, exponential ZNE
brings the error down to about 0.03. Linear ZNE removes only 17–26% of the error.

### 2. Exponential ZNE stays accurate at long times, where linear ZNE is biased

At h/J = 1, the error of linear ZNE grows with time like the raw data, from 0.034 to 0.195.
Exponential ZNE stays between 0.021 and 0.042 through t = 2, where the most amplified
circuit has 2400 two-qubit gates. For example, at t = 2 it gives Mz = 0.517 ± 0.027
against the exact 0.500. Linear ZNE has small error bars but sits 7–45 standard deviations
from the exact result at t ≥ 1.5, and it never lands above the exact value (0 of 36 points).

**Conclusion:** the failure the previous project observed after t ≈ 1.5 is a property of
the linear model, not of the hardware. An extrapolation that follows the exponential
decay of the signal removes it.

### 3. The result is validated and cheap

- **Reproducible:** two runs a day apart both give 0.03. Raw values drift by 0.002–0.016
  between runs, which we report as an extra ±0.01 uncertainty.
- **The model fits:** with five noise levels, the exponential fit passes a χ² test
  (χ²/dof ≤ 2.7) at 19 of 24 points.
- **The residual is hardware noise:** at h/J = 2 the noiseless Trotter circuit is within
  0.011 of exact diagonalization, so discretization is not the limiting factor.
- **Low cost:** each full sweep used 39–58 s of QPU time.

**Conclusion:** the improvement is statistically supported, reproducible and inexpensive.
Our plan uses more quantum resources per experiment than the previous one (12 qubits, five
noise levels, 4000 shots per circuit), but no shot is discarded and each experiment needs
about one minute of real QPU time.

### Limitations

A floor of about 0.02–0.03 remains from noise that CZ folding does not amplify
(single-qubit gates, idling, residual readout error). The exponential model fails its fit
test at 5 of 24 points. The study uses one device, 12 spins and four time points, and no
quantum advantage is claimed: 12 spins are still exactly solvable, which is what lets us
check the results.

## Repository

| Path | Contents |
|---|---|
| `hardware.py` | Backend selection, layout search, folding, Runtime execution ([guide](docs/HARDWARE.md)) |
| `src/qfest/` | TFIM circuits and exact diagonalization |
| `experiments/zne_time_sweep_ibm.py` | Plan, submit and analyze the IBM hardware sweeps |
| `experiments/plot_zne_summary.py` | Summary figures |
| `results/` | Hardware data (JSON) and figures |
| `docs/RESULTS.md` | Full results write-up with references |

Reproduce the analysis from the saved data, without a QPU:

```bash
python experiments/zne_time_sweep_ibm.py analyze --from-json results/zne_time_sweep_ibm_h1_ibm_quebec_lam1-2-3-4-5.json
python experiments/plot_zne_summary.py
```

This project is developed for the University of Ottawa's Qiskit Fall Fest 2026.
