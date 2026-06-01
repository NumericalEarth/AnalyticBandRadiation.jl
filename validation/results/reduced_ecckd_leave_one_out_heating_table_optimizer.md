# Reduced ecCKD Leave-One-Out Heating Table Optimizer

Status: **leave_one_out_heating_table_optimizer_rejected**

| Field | Value |
|---|---:|
| Omitted SW g-point | 25 |
| Candidate count | 16 |
| Probe step | 0.00048828125 |
| Max log scale | 0.001953125 |
| Boundary cap | 0.3 W m^-2 |
| Base objective | 12.332454488 |
| Base boundary forcing | 0.0242811668201 W m^-2 |
| Base heating RMSE | 0.616622724398 K day^-1 |
| Best objective | 12.3659981619 |
| Best objective reduction | -0.0335436739048 |
| Best boundary forcing | 3.69010246526 W m^-2 |
| Best heating RMSE | 0.618299908093 K day^-1 |
| Accepted | false |
| Accepted ridge lambda | n/a |
| Accepted objective | 12.332454488 |
| Accepted boundary forcing | 0.0242811668201 W m^-2 |
| Accepted heating RMSE | 0.616622724398 K day^-1 |
| Accepted moves | 0 |

This diagnostic starts from saved 32x31 leave-one-out weight-refit
models and fits local table-entry moves against the same heating-profile
plus boundary residual used by the retained 16-SW heating optimizer.
It accepts only an exact hard-objective improvement that keeps the
TOA/surface boundary forcing below the hard 0.3 W m^-2 cap.

## Compared Leave-One-Out Rows

| Omitted SW g-point | Status | Base objective | Best objective | Best reduction | Base boundary | Best boundary | Base heating RMSE | Best heating RMSE | Accepted |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | leave_one_out_heating_table_optimizer_rejected | 12.332454488 | 12.3659981619 | -0.0335436739048 | 0.0242811668201 | 3.69010246526 | 0.616622724398 | 0.618299908093 | false |
| 23 | leave_one_out_heating_table_optimizer_rejected | 1.62442314379 | 12.4438908196 | -10.8194676758 | 0.0751392004466 | 3.73316724587 | 0.0812211571893 | 0.0932214516929 | false |

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Boundary forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|---:|
| 1e-10 | 12.3660001975 | 0.57334711501 | 3.69027073923 | 3.69027073923 | 0.618300009874 | false |
| 1e-08 | 12.3660001975 | 0.57334711501 | 3.69027073923 | 3.69027073923 | 0.618300009874 | false |
| 1e-06 | 12.3659995626 | 0.573371217569 | 3.6902878585 | 3.6902878585 | 0.618299978129 | false |
| 0.0001 | 12.3659995848 | 0.573399074389 | 3.69024520611 | 3.69024520611 | 0.618299979241 | false |
| 0.01 | 12.3659990665 | 0.573489492692 | 3.69013584646 | 3.69013584646 | 0.618299953325 | false |
| 1 | 12.3659988948 | 0.573626224592 | 3.6898601393 | 3.6898601393 | 0.618299944739 | false |
| 100 | 12.3659981619 | 0.573538176977 | 3.69010246526 | 3.69010246526 | 0.618299908093 | false |
| 10000 | 12.3659995691 | 0.573394401361 | 3.69024900825 | 3.69024900825 | 0.618299978453 | false |
