# ecCKD Candidate Table Humidity-Split Descent

Status: **humidity_split_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_split_descent_minimax_sw32_candidate.nc`

## Descent Summary

| Metric | Value |
|---|---:|
| Mode | minimax |
| Present scenarios | 3/3 |
| Tested moves | 48 |
| Accepted moves | 4 |
| Initial aggregate loss | 807.1865564430091 |
| Final aggregate loss | 808.4240334667119 |
| Aggregate loss reduction factor | 0.9984692723465973 |
| Initial worst loss ratio | 0.9961868084392239 |
| Final worst loss ratio | 0.9960172397900374 |
| Worst loss ratio reduction factor | 1.000170246700973 |
| Final improved scenarios | 3 |

## Final Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 251.09721910389098 | 0.7634553816868372 | true |
| rel-415 | 330.40047827264766 | 252.56705430779843 | 0.7644270239202837 | true |
| rel-1120 | 1929.2917333373919 | 1921.607826988446 | 0.9960172397900374 | true |

This repeats exact written/reread dry/moist H2O split moves. It is a humidity-aware table optimizer diagnostic, not published-model recovery.
