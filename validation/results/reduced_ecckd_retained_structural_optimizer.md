# Reduced ecCKD Retained Structural Optimizer

Status: **constrained_table_optimizer_improved**

| Field | Value |
|---|---:|
| Base mode | boundary_table_continuation |
| Candidate scope | all_global_residual_probe |
| Residual mode | all_shortwave |
| Include Rayleigh candidates | true |
| Candidate count | 12 |
| Probe step | 0.001953125 |
| Max log scale | 0.03125 |
| Base objective | 7.40057407226 |
| Best ridge lambda | 0.0001 |
| Best exact objective | 7.25884945943 |
| Best objective reduction | 0.141724612829 |
| Best TOA forcing error | 2.17765483783 W m^-2 |
| Best surface forcing error | 1.98402949915 W m^-2 |
| Accepted moves | 12 |
| Accepted | true |

This retained optimizer starts from the current boundary-table-continuation
reduced model and records the accepted all-shortwave residual constrained
table moves. The canonical reduced accuracy path consumes these moves after
the pressure-temperature retained chain.
