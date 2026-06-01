# ecCKD Candidate Table Written Minimax Descent

Status: **written_minimax_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_written_minimax_sw32_candidate.nc`

## Written Minimax Descent

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Tested moves | 144 |
| Accepted moves | 6 |
| Initial aggregate loss | 806.8325371593586 |
| Final aggregate loss | 848.1447379487423 |
| Aggregate loss reduction factor | 0.9512910957989338 |
| Initial worst loss ratio | 0.9962357775961005 |
| Final worst loss ratio | 0.9934182526784027 |
| Worst loss ratio reduction factor | 1.0028361920169087 |
| Final improved scenarios | 3 |

## Final Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 313.2035809146191 | 0.9522883617198508 | true |
| rel-415 | 330.40047827264766 | 314.6370102926893 | 0.9522898148865563 | true |
| rel-1120 | 1929.2917333373919 | 1916.5936226389185 | 0.9934182526784027 | true |

This repeats exact written-file coordinate scans, accepting moves by the reread worst scenario loss ratio. It targets the high-humidity representative state and is not published-model recovery.
