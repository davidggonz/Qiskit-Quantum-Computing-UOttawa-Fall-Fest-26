# Results log

All numbers use the previous report's error convention: 100 * max_t |O - O_exact| / max_t |O_exact|.
Simulations use FakeMarrakesh noise unless stated. Scripts in `experiments/`, data in `results/`.

## 1. Table 2 reproduced exactly (Trotter vs ED, noiseless)
`experiments/table2_trotter.py` matches every entry of the old report's Table 2. See also `docs/TABLE1_DISCREPANCIES.md`.

## 2. MPF vs Trotter at equal max depth, shot noise only
`experiments/mpf_short_times.py`. N=6, k_max=8, t <= 2, 10k shots: best MPF set halves the error at h/J <= 1
([4,6,8] with ||x||_1 <= 1.5) and gets h/J=2 under 5% ([6,8] exact). At 1k shots plain Trotter wins.

## 3. Zero-SWAP layout and better ZNE fits (Boutaina's canonical run vs the ported pipeline)
`experiments/compare_hardware_runs.py`, figures `results/figs/compare_*.png`. Signed <Z>, h/J=0.5, 20 steps:
N=12 ring with the same 480 CZ as her N=6 ring (with SWAPs), half the depth; raw 0.699 vs 0.449; error after ZNE 1.7%
(exponential) vs 40% (her linear). Exact 0.868 for both.

## 4. MPF + ZNE under noise (N=6 open chain, k_max=8, t in (0,2], 8 time points)
`experiments/mpf_zne_noisy.py`, figure `results/figs/mpf_zne_noisy_n6_chain.png`. Mean +- std over seeds.

| h/J | obs | shots/circuit | linear ZNE (old) | Richardson ZNE | exp ZNE | MPF + exp ZNE | MPF + Rich ZNE | ideal Trotter |
|---|---|---|---|---|---|---|---|---|
| 0.5 | Mzz | 8k (5 seeds) | 4.27 ± 0.40 | 3.20 ± 0.48 | **3.09 ± 0.45** | 3.74 ± 0.67 | 3.75 ± 0.80 | 2.27 |
| 2 | Mzz | 8k (5 seeds) | 4.93 ± 0.32 | **3.48 ± 0.50** | 3.77 ± 0.40 | 6.16 ± 1.65 | 7.07 ± 1.92 | 4.17 |
| 2 | Mzz | 128k (3 seeds) | 5.02 ± 0.20 | 3.05 ± 0.26 | **3.04 ± 0.12** | 3.91 ± 0.45 | 4.06 ± 0.41 | 4.17 |
| 2 | Zavg | 8k (5 seeds) | **8.37 ± 1.03** | 9.10 ± 2.08 | 10.27 ± 2.95 | 11.17 ± 5.68 | 10.70 ± 2.92 | 8.55 |
| 2 | Zavg | 32k (3 seeds) | 8.39 ± 0.32 | 9.21 ± 0.83 | 9.50 ± 1.04 | **4.95 ± 1.29** | 5.94 ± 1.50 | 8.55 |
| 2 | Zavg | 128k (3 seeds) | 8.27 ± 0.08 | 8.43 ± 0.48 | 8.57 ± 0.50 | 6.12 ± 2.95 | **3.65 ± 0.46** | 8.55 |

Conclusions:
- **Exponential / Richardson ZNE beat the old linear ZNE for Mzz in every setting** (about 1.3-1.6x lower error).
- **ZNE alone cannot remove Trotter error**: it recovers the ideal Trotter value, so where Trotter error dominates
  (Zavg at h/J=2, 8.6%) all Trotter+ZNE variants stay around 8-10%.
- **MPF + ZNE is shot-noise limited**: ZNE and MPF both amplify statistical noise (MPF by ||x||_1 = 3.6 at h/J=2).
  At 8k shots/circuit it loses; with >= 32k shots/circuit it is the only method that removes the Trotter error
  (Zavg h/J=2: 8.4% -> 3.7-5.0%). For Mzz, where noise dominates, it does not pay off.
- The exponential fit needed a guard: with near-zero data it diverged (600%-15000%); it now falls back to
  Richardson/linear when |estimate| > 1 (`qfest.extrapolation.exponential`).

Implication for hardware: use Trotter + exponential/Richardson ZNE as the main method. Run MPF + ZNE only for one
showcase case where Trotter error dominates (Zavg, h/J=2), with >= 32k shots per circuit and few time points.
