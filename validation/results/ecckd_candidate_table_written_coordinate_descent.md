# ecCKD Candidate Table Written Coordinate Descent

Status: **written_coordinate_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_written_coordinate_descent_sw32_candidate.nc`

## Written Coordinate Descent

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Iterations | 6 |
| Accepted moves | 6 |
| Initial aggregate loss | 808.690677162573 |
| Final aggregate loss | 806.8325371593586 |
| Aggregate loss reduction factor | 1.0023030057882352 |
| Initial worst loss ratio | 0.9959808898432916 |
| Final worst loss ratio | 0.9962357775961005 |
| Final improved scenarios | 3 |

## Final Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 248.5080485320662 | 0.7555830674723697 | true |
| rel-415 | 330.40047827264766 | 249.96011277490453 | 0.7565367764650648 | true |
| rel-1120 | 1929.2917333373919 | 1922.029450171105 | 0.9962357775961005 | true |

This repeats exact written-file coordinate scans, accepting a move only when the reread multi-sample aggregate loss improves. It is a written-file optimizer diagnostic, not published-model recovery.
