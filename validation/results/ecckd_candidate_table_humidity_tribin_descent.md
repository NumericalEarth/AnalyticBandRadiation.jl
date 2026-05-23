# ecCKD Candidate Table Humidity-Tribin Descent

Status: **humidity_tribin_descent_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_humidity_tribin_descent_sw32_candidate.nc`

## Descent Summary

| Metric | Value |
|---|---:|
| Mode | aggregate |
| Present scenarios | 3/3 |
| Tested moves | 72 |
| Accepted moves | 4 |
| Initial aggregate loss | 804.9001781269905 |
| Final aggregate loss | 803.5814834273937 |
| Aggregate loss reduction factor | 1.0016410217592027 |
| Initial worst loss ratio | 0.9965216904173317 |
| Final worst loss ratio | 0.9967387288858741 |
| Worst loss ratio reduction factor | 0.9997822513942195 |
| Final improved scenarios | 3 |

## Final Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 243.16507933329618 | 0.7393378911874813 | true |
| rel-415 | 330.40047827264766 | 244.57958101214786 | 0.740251897608695 | true |
| rel-1120 | 1929.2917333373919 | 1922.9997899367368 | 0.9967387288858741 | true |

This repeats exact written/reread low/mid/high H2O split moves. It is a structured humidity-axis table optimizer diagnostic, not published-model recovery.
