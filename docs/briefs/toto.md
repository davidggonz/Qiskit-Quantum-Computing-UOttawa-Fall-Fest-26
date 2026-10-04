# Toto: noise models, plots and docs

**You own:** `src/qfest/noise.py`, `src/qfest/plots.py`, the README and the slides.

## The idea in plain words
- Real quantum chips are noisy. Gates sometimes apply a small random error, and measurements sometimes read 0 as 1.
- A **noise model** reproduces that noise in a simulator (Qiskit Aer), so we can test ideas without spending scarce hardware time.
- Your plots are what the judges see. They must make the comparison with the old project obvious.

## Steps
1. **`simple_model()` (day 1).** Use `qiskit_aer.noise`: `depolarizing_error` for 1- and 2-qubit gates plus a `ReadoutError`. Start with p1q = 4e-4, p2q = 3e-3, p_meas = 3e-3.
   *Test:* a 2-qubit circuit run with noise gives slightly wrong counts, and without noise gives exact ones.
2. **`from_backend()` (day 1–2).** Build a noise model from Boutaina's fake IBM backend (`AerSimulator.from_backend`).
3. **Plot functions (day 2).** Write one function per figure, each reading `results/*.json` (schema in `results.py`). Make a small fake JSON to develop against. The figures are:
   - observable vs time: ED vs Trotter vs MPF vs mitigated
   - error vs time on a log scale
   - Iceberg discard rate vs circuit depth
   - SWAP overhead, all-to-all vs heavy-hex (bar chart)
4. **Summary table (day 3).** Show the old report's numbers next to ours, for example "8.46% → ?" for Trotter error at h/J = 2 and "42.3% → ?" for ZNE.
5. **README and slides (day 3–4).** Write how to reproduce our results. The slide story is problem → old limitations → our fixes → results → what's next.

Use clear axis labels and colour-blind-friendly colours, and save the figures to `results/figs/`.
