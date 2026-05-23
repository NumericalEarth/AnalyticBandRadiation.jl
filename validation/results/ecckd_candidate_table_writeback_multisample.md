# ecCKD Candidate Table-Writeback Multisample

Status: **table_writeback_multisample_passed**

## Blockers

None.

## Scenario Scores

| Scenario | Present | Baseline loss | Candidate loss | Loss ratio | Improved |
|---|---:|---:|---:|---:|---:|
| rel-180 | true | 328.8957352675901 | 252.03195214968295 | 0.7662974162453926 | true |
| rel-415 | true | 330.40047827264766 | 253.5083950011559 | 0.76727611390429 | true |
| rel-1120 | true | 1929.2917333373919 | 1921.6807505232734 | 0.9960550378760228 | true |

This rescans the written table-continuation candidate across multiple representative CKDMIP SW scenarios. It checks generalization of the writeback update; it is not full published-model recovery.
