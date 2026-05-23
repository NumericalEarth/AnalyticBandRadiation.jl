# ecCKD Candidate Table Humidity-Tribin Probe

Status: **humidity_tribin_probe_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_tribin_sw32_candidate.nc`

## Probe Summary

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Tested moves | 18 |
| Seed aggregate loss | 805.1956684744637 |
| Best aggregate loss | 804.9001781269905 |
| Aggregate loss reduction factor | 1.0003671142776498 |
| Seed worst loss ratio | 0.9964754504419561 |
| Best minimax worst loss ratio | 0.9964319180920038 |
| Worst loss ratio reduction factor | 1.0000436882331465 |

This writes and rereads one-step low/mid/high H2O split-table moves. It tests whether a three-bin humidity structure decouples representative states better than the dry/moist split; it is not published-model recovery.
