# Reduced ecCKD Leave-One-Out Weight Coordinate Descent Continuation

Status: **weight_coordinate_descent_continuation_improved**

| Field | Value |
|---|---:|
| Omitted SW g-point | 23 |
| Model | 32x31 |
| Coordinate count | 31 |
| Max iterations | 8 |
| Completed iterations | 8 |
| Accepted moves | 8 |
| Initial objective | 1.26634555306 |
| Final objective | 1.08338530362 |
| Objective reduction | 0.18296024944 |
| Final TOA forcing | 0.299433437452 W m^-2 |
| Final surface forcing | 0.278654452263 W m^-2 |
| Final boundary forcing | 0.299433437452 W m^-2 |
| Final heating RMSE | 0.054169265181 K day^-1 |

This diagnostic continues from the saved six-step exact
weight-coordinate descent using a smaller log-scale ladder and all
31 remaining shortwave weights. It tests whether the boundary-tight
32x31 model still has safe quadrature-weight descent directions.

## Iterations

| Iteration | Accepted | G-point | Log scale | Objective | Reduction | Boundary forcing | Heating RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | true | 25 | 0.00316227766017 | 1.24510646479 | 0.0212390882714 | 0.299946497208 | 0.0622553232394 |
| 2 | true | 25 | 0.001 | 1.23856148974 | 0.00654497504545 | 0.299969320306 | 0.0619280744871 |
| 3 | true | 25 | 0.001 | 1.23210179497 | 0.00645969477393 | 0.299992166238 | 0.0616050897484 |
| 4 | true | 25 | 0.000316227766017 | 1.23007708458 | 0.00202471039252 | 0.299999395511 | 0.0615038542288 |
| 5 | true | 16 | -0.001 | 1.22867698481 | 0.00140009976822 | 0.290898628449 | 0.0614338492404 |
| 6 | true | 32 | 0.00316227766017 | 1.12529293058 | 0.103384054228 | 0.296880867243 | 0.056264646529 |
| 7 | true | 32 | 0.001 | 1.10304433733 | 0.0222485932477 | 0.299361018706 | 0.0551522168666 |
| 8 | true | 25 | 0.00316227766017 | 1.08338530362 | 0.0196590337126 | 0.299433437452 | 0.054169265181 |
