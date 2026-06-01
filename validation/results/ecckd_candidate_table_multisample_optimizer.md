# ecCKD Candidate Table Multisample Optimizer

Status: **table_multisample_optimizer_passed**

## Blockers

None.

Candidate: `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_table_multisample_sw32_candidate.nc`

## Aggregate Optimizer

| Metric | Value |
|---|---:|
| Present scenarios | 3/3 |
| Iterations | 2 |
| Accepted steps | 2 |
| Initial aggregate loss | 809.0736990997606 |
| Final aggregate loss | 808.9416235913606 |
| Aggregate loss reduction factor | 1.0001632695172906 |
| Best written checkpoint | 2 |
| Selected written aggregate loss | 808.9416241974774 |
| Writeback status | table_writeback_multisample_passed |
| Writeback improved scenarios | 3 |
| Writeback worst loss ratio | 0.9959467737433978 |

## Written Candidate Scenario Scores

| Scenario | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|
| rel-180 | 328.8957352675901 | 251.93871537566943 | 0.7660139319553404 | true |
| rel-415 | 330.40047827264766 | 253.41427978957984 | 0.7669912619813506 | true |
| rel-1120 | 1929.2917333373919 | 1921.4718774271832 | 0.9959467737433978 | true |

This optimizes four real ecCKD shortwave table multipliers against an aggregate projected original-objective loss over multiple CKDMIP SW scenarios, writes a CKD-definition candidate, and rescans the written file. It is still a small-parameter optimizer probe, not published-model recovery.
