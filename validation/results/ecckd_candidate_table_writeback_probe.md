# ecCKD Candidate Table-Writeback Probe

Status: **table_writeback_probe_passed**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_scaled_sw32_candidate.nc`

## Writeback Probe

| Metric | Value |
|---|---:|
| Baseline SW projected loss | 330.40047827264766 |
| Written SW projected loss | 264.1244985078305 |
| Loss reduction factor | 1.2509270444023286 |
| Accepted writeback | true |
| Baseline non-finite flux count | 0 |
| Written non-finite flux count | 0 |
| Candidate metrics status | failed |
| Candidate worst log-coefficient RMSE | 0.0038167317109358582 |

This writes the optimized table multipliers into a CKD-definition NetCDF candidate and rescans that file through the normal reader and projected original-objective score. It is a writeback probe, not full published-model recovery.
