# Reduced ecCKD Retained Boundary-Probe Expansion

Status: **boundary_probe_expansion_rejected**

| Field | Value |
|---|---:|
| Configurations | 3 |
| Best label | boundary_probe_32_small_tiny |
| Best exact objective | 7.19450656389 |
| Best objective reduction | -0.000982379344379 |
| Best TOA forcing | 2.15835196917 W m^-2 |
| Best surface forcing | 2.02981871317 W m^-2 |
| Any accepted | false |

This diagnostic uses the current retained base with a TOA+surface
boundary residual and accepts only strict Pareto-safe full-objective
updates that do not regress either TOA or surface forcing.

## Configurations

| Label | Scope | Candidates | Probe step | Max log scale | Base objective | Base TOA | Base surface | Best objective | Reduction | TOA | Surface | Pareto safe | Accepted |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| boundary_probe_32_medium_tiny | all_global_objective_probe | 32 | 0.000244140625 | 0.001953125 | 7.19352418454 | 2.15805725536 | 2.02792974588 | 7.19450661621 | -0.000982431662588 | 2.15835198486 | 2.02981881105 | false | false |
| boundary_probe_32_small_tiny | all_global_objective_probe | 32 | 0.0001220703125 | 0.0009765625 | 7.19352418454 | 2.15805725536 | 2.02792974588 | 7.19450656389 | -0.000982379344379 | 2.15835196917 | 2.02981871317 | false | false |
| boundary_probe_48_residual_small | all_global_residual_probe | 48 | 0.0009765625 | 0.0078125 | 7.19352418454 | 2.15805725536 | 2.02792974588 | 7.194506929 | -0.000982744453685 | 2.1583520787 | 2.02981905282 | false | false |
