# David: team leader, physics lead, integrator

**You own:** `src/qfest/tfim.py`, `src/qfest/mpf.py`, code review, the final analysis and write-up.

## Steps
1. **Kick-off (day 1, 30 min).** Share the briefs, confirm everyone ran `pytest`, and confirm the hackathon's IBM QPU allocation.
2. **Validate Trotter against the old report.** For N = 6 and T = 20, compute the max Mzz deviation versus ED for Δt in {0.2, 0.1, 0.05, 0.025} and h/J in {0.5, 1, 2}. It should match the report's Table 2 (4.47% … 0.073% at h = 0.5).
   *Done when:* the table is reproduced and saved with `qfest.results.save`.
3. **Multi-product formulas** in `mpf.py`. Start with Richardson coefficients for step counts like [1, 2, 4], then compare with `qiskit-addon-mpf`.
   *Done when:* MPF beats plain 2nd-order Trotter at equal max depth for h/J = 2. The target is below 5%; the report had 8.46%.
4. **Freeze the hardware circuits** (end of day 2) with Boutaina. Keep it small: one N, one h, a few times, plus ZNE noise factors.
5. **Day 3 hardware run** with Boutaina. You shadow her and take over if she runs out of hours.
6. **Review PRs** daily. Check that tests exist and the results JSON follows the schema in `results.py`.
7. **Write-up (day 4).** Compare against the old report: Trotter error, ZNE vs better extrapolators, discard rates, circuit cost. Mention the corrected Table 1 (`docs/TABLE1_DISCREPANCIES.md`).

**Stretch:** 2×2 Hubbard in simulation only.
