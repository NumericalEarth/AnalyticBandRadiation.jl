# ecCKD Candidate Table Humidity-Tribin Weighted Descent

Status: **humidity_tribin_weighted_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_tribin_weighted_descent_sw32_candidate.nc`

## Best Weighted Descent Summary

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Tested weight sets | 4 |
| Best scenario weights | Dict("rel-1120" => 2.0) |
| Best tested moves | 72 |
| Best accepted moves | 4 |
| Best weighted loss ratio | 0.8682668116419812 |
| Weighted loss ratio reduction factor | 1.0036961710911658 |
| Best aggregate loss | 803.5814834273937 |
| Aggregate loss reduction factor | 1.0016410217592027 |
| Best worst loss ratio | 0.9967387288858741 |
| Worst loss ratio reduction factor | 0.9997822513942195 |
| Final improved scenarios | 3 |

## Best Final Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 243.16507933329618 | 0.7393378911874813 | true |
| rel-415 | 330.40047827264766 | 244.57958101214786 | 0.740251897608695 | true |
| rel-1120 | 1929.2917333373919 | 1922.9997899367368 | 0.9967387288858741 | true |

## Weight-Set Sweep

| Scenario weights | Accepted moves | Weighted ratio | Worst ratio | Aggregate loss |
|---|---:|---:|---:|---:|
| Dict("rel-1120" => 2.0) | 4 | 0.8682668116419812 | 0.9967387288858741 | 803.5814834273937 |
| Dict("rel-1120" => 4.0) | 4 | 0.9110907840566121 | 0.9967387288858741 | 803.5814834273937 |
| Dict("rel-1120" => 8.0) | 4 | 0.945349961988317 | 0.9967387288858741 | 803.5814834273937 |
| Dict("rel-1120" => 16.0) | 4 | 0.9681894139427869 | 0.9967387288858741 | 803.5814834273937 |

This exact written/reread scan optimizes a weighted mean of per-scenario loss ratios, emphasizing high-humidity rel-1120 without using raw-loss scale as the only objective. It is a diagnostic for objective weighting, not published-model recovery.
