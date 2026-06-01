# ecCKD Candidate Table Humidity-Split Descent

Status: **humidity_split_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_split_descent_sw32_candidate.nc`

## Descent Summary

| Metric | Value |
|---|---:|
| Mode | aggregate |
| Present scenarios | 3/3 |
| Tested moves | 48 |
| Accepted moves | 4 |
| Initial aggregate loss | 806.460688875991 |
| Final aggregate loss | 804.8071816437472 |
| Aggregate loss reduction factor | 1.002054538366403 |
| Initial worst loss ratio | 0.996287642294726 |
| Final worst loss ratio | 0.9965278110668085 |
| Worst loss ratio reduction factor | 0.9997589944109784 |
| Final improved scenarios | 3 |

## Final Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 245.1998108665163 | 0.7455244461197447 | true |
| rel-415 | 330.40047827264766 | 246.628866132725 | 0.7464543254359516 | true |
| rel-1120 | 1929.2917333373919 | 1922.592867932 | 0.9965278110668085 | true |

This repeats exact written/reread dry/moist H2O split moves. It is a humidity-aware table optimizer diagnostic, not published-model recovery.
