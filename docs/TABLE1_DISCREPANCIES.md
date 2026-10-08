# Discrepancies in Table 1 of the previous report

Our ED (`src/qfest/ed.py`, periodic chain, J = 1, same observable definitions) disagrees with
several entries of the report's Table 1. Our values satisfy exact identities that the report's
values violate: Kramers-Wannier duality `Mx(h) = Mzz(1/h)` and self-duality `Mx = Mzz` at h = 1.

| N | h/J | quantity | report | ours |
|---|---|---|---|---|
| 6 | 0.5 | Mx  | 0.251951 | 0.265198 (= report's own Mzz at h=2) |
| 6 | 1.0 | Mx  | 0.643915 | 0.643951 |
| 6 | 1.0 | Mzz | 0.645951 | 0.643951 |
| 6 | 2.0 | Mz  | 0.548173 | 0.554173 |
| 8 | 0.5 | Mz  | 0.968903 | 0.969303 |
| 8 | 1.0 | Mz  | 0.875827 | 0.785827 |
| 8 | 1.0 | Mzz | 0.406729 | 0.640729 |
| 10 | 0.5 | Mz | 0.967859 | 0.968589 |

These look like transcription errors in the report's table (digits transposed or mistyped), not
different physics. All other entries match to 6 decimals. Use our ED values as the reference and
mention the correction in the final write-up.
