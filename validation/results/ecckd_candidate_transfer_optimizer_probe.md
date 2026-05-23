# ecCKD Candidate Transfer Optimizer Probe

Status: **optimizer_probe_passed**

## Blockers

None.

## Optimizer Probe

| Metric | Value |
|---|---:|
| Parameters | 4 |
| Initial projected loss | 686.9810686483928 |
| Final projected loss | 620.0785606939874 |
| Loss reduction factor | 1.107893599610231 |
| Accepted step | true |
| Step size | 0.001 |
| Gradient method | central_finite_difference_on_band_projected_transfer_loss |
| Gradient norm | 314.73818946962643 |

This optimizes four global log-scale transfer parameters against the candidate-driven, g-point-to-CKDMIP-band projected original-objective loss. It is an optimizer probe for the transfer path, not a published CKD-definition recovery.
