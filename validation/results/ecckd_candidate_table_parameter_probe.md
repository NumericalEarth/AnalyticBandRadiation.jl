# ecCKD Candidate Table-Parameter Probe

Status: **table_parameter_probe_passed**

## Blockers

None.

## Optimizer Probe

| Metric | Value |
|---|---:|
| Parameters | 4 |
| Initial SW projected loss | 330.40047827264806 |
| Final SW projected loss | 264.1244998704228 |
| Loss reduction factor | 1.2509270379489206 |
| Accepted step | true |
| Step size | 0.001 |
| Gradient method | central_finite_difference_on_in_memory_ckd_table_scales |
| Gradient norm | 313.7507279595536 |

This optimizes four real in-memory ecCKD shortwave table multipliers against the candidate-driven projected CKDMIP-band loss. It is a table-parameter optimizer probe, not a full CKD-definition writeback or published-model recovery.
