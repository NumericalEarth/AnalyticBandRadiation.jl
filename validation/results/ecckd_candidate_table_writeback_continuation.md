# ecCKD Candidate Table-Writeback Continuation

Status: **table_writeback_continuation_passed**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_scaled_sw32_continuation_candidate.nc`

## Continuation

| Metric | Value |
|---|---:|
| Iterations | 4 |
| Accepted steps | 4 |
| Initial in-memory loss | 330.40047827264806 |
| Final in-memory loss | 253.5083948490462 |
| In-memory loss reduction factor | 1.3033117836961097 |
| Baseline written SW projected loss | 330.40047827264766 |
| Continuation written SW projected loss | 253.5083950011559 |
| Writeback loss reduction factor | 1.303311782914097 |
| Best written checkpoint | 4 |
| Accepted writeback | true |
| Candidate metrics status | failed |
| Candidate worst log-coefficient RMSE | 0.0151412979344743 |

This runs a short continuation on real in-memory ecCKD table multipliers, writes the final multipliers into a CKD-definition NetCDF candidate, and verifies that the written file still reduces the projected original-objective loss. It is not full published-model recovery.
