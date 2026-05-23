# ecCKD Candidate Table Written Coordinate Scan

Status: **written_coordinate_scan_improved**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_written_coordinate_sw32_candidate.nc`

## Written Coordinate Scan

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Tested moves | 24 |
| Seed aggregate loss | 808.9416241974774 |
| Best aggregate loss | 808.690677162573 |
| Aggregate loss reduction factor | 1.0003103127586246 |
| Seed worst loss ratio | 0.9959467737433978 |
| Best worst loss ratio | 0.9959808898432916 |
| Best improved scenarios | 3 |

## Best Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 251.53077223993682 | 0.764773590132724 | true |
| rel-415 | 330.40047827264766 | 253.00356191110035 | 0.7657481709282542 | true |
| rel-1120 | 1929.2917333373919 | 1921.537697336682 | 0.9959808898432916 | true |

This exact scan writes and rereads one-coordinate table-parameter moves around the selected multi-sample candidate. A move is useful only if the written file improves the multi-sample score; in-memory optimizer loss is not used for acceptance.
