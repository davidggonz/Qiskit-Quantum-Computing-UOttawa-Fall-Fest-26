# Ririsha: classical baseline and numbers

**You own:** `src/qfest/ed.py`, `src/qfest/extrapolation.py`, `src/qfest/metrics.py`.
You can start without knowing Qiskit; these modules are plain Python with numpy and scipy.

## The idea in plain words
- Our physical system (a ring of spins) is described by a big matrix called the **Hamiltonian**. For 12 spins it is 4096 × 4096, which is small enough to solve exactly on a laptop.
- Solving it exactly (**ED**, exact diagonalization) gives the **answer key**. Every quantum result is graded against your numbers.
- **ZNE** (zero-noise extrapolation): we run a circuit at noise levels 1×, 3× and 5×, then fit a curve and read off the value at 0× noise. The old project used a straight line and overshot. You'll provide better curve fits.

## Steps
1. **Setup.** Run `pytest`, then read `ed.py` and `tests/test_ed_and_tfim.py` to see how the answer key is checked.
2. **Make reference data (day 1).** Run `ed.quench` for N = 6, 8, 10, 12 and h = 0.5, 1, 2, with times 0 to 20 in steps of 0.05. Save each with `qfest.results.save` (method `"ed"`).
   *Done when:* the JSON files are in `results/` and Toto can plot them.
3. **Extrapolators (day 2).** Write tests for `linear`, `richardson` and `exponential`:
   - Fake data where you know the true answer.
   - Data with small random noise added.
   - A case where the exponential fit fails, so the code must fall back to linear.

   Fix any bugs you find.
4. **Metrics.** Add tests for `metrics.py`, and add a function that returns error vs time for a whole curve.
5. **Compare extrapolators (day 3).** On the noisy-simulation data from Eli and Toto, compute which extrapolator gets closest to ED and save the comparison.
6. **If you have time:** IBM Quantum Learning, "Basics of quantum information", lesson 1.

Ask David about physics questions and Toto about plot format.
