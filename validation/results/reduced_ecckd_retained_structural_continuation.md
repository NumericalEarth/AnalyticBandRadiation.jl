# Reduced ecCKD Retained Structural Continuation

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_continuation |
| Candidate scope | all_global_residual_probe |
| Residual mode | all_shortwave |
| Include Rayleigh candidates | true |
| Candidate count | 12 |
| Probe step | 0.0009765625 |
| Max log scale | 0.0078125 |
| Base objective | 7.25884945943 |
| Best ridge lambda | 0.0001 |
| Best exact objective | 7.22363619225 |
| Best objective reduction | 0.0352132671795 |
| Best TOA forcing error | 2.16709085767 W m^-2 |
| Best surface forcing error | 1.975674676 W m^-2 |
| Accepted moves | 12 |
| Accepted | true |

This retained continuation starts from the current reduced model after
the retained structural optimizer and records a smaller-trust-region
all-shortwave residual constrained table continuation.
