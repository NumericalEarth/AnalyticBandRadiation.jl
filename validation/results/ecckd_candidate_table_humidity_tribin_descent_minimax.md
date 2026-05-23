# ecCKD Candidate Table Humidity-Tribin Descent

Status: **humidity_tribin_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_tribin_descent_minimax_sw32_candidate.nc`

## Descent Summary

| Metric | Value |
|---|---:|
| Mode | minimax |
| Present scenarios | 3/3 |
| Tested moves | 72 |
| Accepted moves | 4 |
| Initial aggregate loss | 805.4766666231361 |
| Final aggregate loss | 806.453319552341 |
| Aggregate loss reduction factor | 0.9987889529306582 |
| Initial worst loss ratio | 0.9964319180920038 |
| Final worst loss ratio | 0.9962810629818464 |
| Worst loss ratio reduction factor | 1.0001514182249995 |
| Final improved scenarios | 3 |

## Final Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 247.89762003519377 | 0.7537270735161824 | true |
| rel-415 | 330.40047827264766 | 249.34551973036372 | 0.7546766307178372 | true |
| rel-1120 | 1929.2917333373919 | 1922.1168188914658 | 0.9962810629818464 | true |

This repeats exact written/reread low/mid/high H2O split moves. It is a structured humidity-axis table optimizer diagnostic, not published-model recovery.
