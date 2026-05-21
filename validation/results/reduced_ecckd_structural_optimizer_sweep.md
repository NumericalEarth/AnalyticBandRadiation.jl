# Reduced ecCKD Structural Optimizer Sweep

Status: **structural_optimizer_sweep_improved**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | retained_all_shortwave_residual_probe |
| Best base objective | 7.40057407226 |
| Best exact objective | 7.25884945943 |
| Best objective reduction | 0.141724612829 |
| Best TOA forcing | 2.17765483783 W m^-2 |
| Best surface forcing | 1.98402949915 W m^-2 |

This diagnostic runs the existing constrained table optimizer against
the retained boundary-table-continuation model with multiple structural
candidate and residual definitions. It does not overwrite the retained
optimizer artifact.

## Configurations

| Label | Scope | Residual | Base objective | Best objective | Reduction | TOA | Surface | Accepted |
|---|---|---|---:|---:|---:|---:|---:|---:|
| retained_all_shortwave_residual_probe | all_global_residual_probe | all_shortwave | 7.40057407226 | 7.25884945943 | 0.141724612829 | 2.17765483783 | 1.98402949915 | true |
| retained_boundary_residual_probe | all_global_residual_probe | boundary | 7.40057407226 | 7.60073054281 | -0.200156470552 | 2.28021916284 | 2.0915892285 | false |
