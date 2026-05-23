# ecCKD Candidate Table Humidity-Tribin Constrained Descent

Status: **humidity_tribin_constrained_descent_no_descent**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_tribin_constrained_descent_sw32_candidate.nc`

## Best Constrained Descent Summary

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Tested tolerances | 4 |
| Best worst-ratio tolerance | 0.0 |
| Best tested moves | 18 |
| Best accepted moves | 0 |
| Best aggregate loss | 804.9001781269905 |
| Aggregate loss reduction factor | 1.0 |
| Best worst loss ratio | 0.9965216904173317 |
| Worst loss ratio reduction factor | 1.0 |
| Final improved scenarios | 3 |

## Best Final Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 245.34475500482273 | 0.7459651454744763 | true |
| rel-415 | 330.40047827264766 | 246.77471996258706 | 0.7468957710132238 | true |
| rel-1120 | 1929.2917333373919 | 1922.5810594135617 | 0.9965216904173317 | true |

## Tolerance Sweep

| Worst-ratio tolerance | Accepted moves | Aggregate loss | Worst ratio |
|---|---:|---:|---:|
| 0.0 | 0 | 804.9001781269905 | 0.9965216904173317 |
| 1.0e-5 | 4 | 804.7212291217243 | 0.9965461799195802 |
| 5.0e-5 | 4 | 804.1543387352618 | 0.9966327936121814 |
| 0.0001 | 4 | 803.5814834273937 | 0.9967387288858741 |

This exact written/reread scan accepts aggregate-loss moves only when the worst scenario loss ratio stays below the incumbent plus a tolerance. It tests a constrained multi-objective rule; it is not published-model recovery.
