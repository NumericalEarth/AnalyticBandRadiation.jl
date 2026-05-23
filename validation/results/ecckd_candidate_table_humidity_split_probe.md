# ecCKD Candidate Table Humidity-Split Probe

Status: **humidity_split_probe_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_split_sw32_candidate.nc`

## Probe Summary

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Tested moves | 12 |
| Seed aggregate loss | 806.8325371593586 |
| Best aggregate loss | 806.460688875991 |
| Aggregate loss reduction factor | 1.0004610866822112 |
| Seed worst loss ratio | 0.9962357775961005 |
| Best minimax worst loss ratio | 0.9961868084392239 |
| Worst loss ratio reduction factor | 1.0000491566004104 |

This writes and rereads one-step dry/moist H2O split-table moves. It tests whether humidity-aware table parameters can decouple low-/mid-humidity and high-humidity representative states; it is not published-model recovery.
