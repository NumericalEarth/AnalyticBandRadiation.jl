# Reduced ecCKD Retained Structural Continuation 4

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_continuation |
| Candidate scope | all_global_residual_probe |
| Residual mode | all_shortwave |
| Include Rayleigh candidates | true |
| Candidate count | 12 |
| Probe step | 0.0001220703125 |
| Max log scale | 0.0009765625 |
| Base objective | 7.19781262031 |
| Best ridge lambda | 0.0001 |
| Best exact objective | 7.19352418454 |
| Best objective reduction | 0.00428843576586 |
| Best TOA forcing error | 2.15805725536 W m^-2 |
| Best surface forcing error | 2.02792974588 W m^-2 |
| Accepted moves | 12 |
| Accepted | true |

This retained continuation starts from the current reduced model after
the third retained structural continuation and records another
smaller-trust-region all-shortwave residual constrained table step.
