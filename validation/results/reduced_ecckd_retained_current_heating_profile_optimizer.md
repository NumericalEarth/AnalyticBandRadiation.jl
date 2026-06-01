# Reduced ecCKD Retained Current Heating-Profile Optimizer

Status: **current_heating_profile_optimizer_rejected**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Candidate count | 8 |
| Probe step | 0.00048828125 |
| Max log scale | 0.001953125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.97212754254 |
| Base TOA forcing | 2.09163826276 W m^-2 |
| Base surface forcing | 2.00694307171 W m^-2 |
| Base heating RMSE | 0.348575997036 K day^-1 |
| Best objective | 9.89161448599 |
| Best objective reduction | -2.91948694345 |
| Best TOA forcing | 2.37987968063 W m^-2 |
| Best surface forcing | 2.9674843458 W m^-2 |
| Best heating RMSE | 0.34485215703 K day^-1 |
| Accepted | false |
| Accepted ridge lambda | n/a |
| Accepted objective | 6.97212754254 |
| Accepted TOA forcing | 2.09163826276 W m^-2 |
| Accepted surface forcing | 2.00694307171 W m^-2 |
| Accepted heating RMSE | 0.348575997036 K day^-1 |
| Accepted moves | 0 |

This diagnostic fits local active table-entry moves against a residual
that combines full heating-rate profile errors with TOA and surface
net-boundary errors. It exact-evaluates every ridge row and accepts
only a hard-objective improvement with non-regressing TOA, capped
surface forcing, and non-regressing heating-rate RMSE.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 9.89278899427 | 2.37942251673 | 2.96783669828 | 0.344851781034 | false |
| 1e-06 | 9.89278899427 | 2.37942251673 | 2.96783669828 | 0.344851781034 | false |
| 0.0001 | 9.89273020415 | 2.37944589265 | 2.96781906124 | 0.344851806513 | false |
| 0.01 | 9.89284799692 | 2.37942437719 | 2.96785439908 | 0.34485185946 | false |
| 1 | 9.89179256092 | 2.3798515067 | 2.96753776828 | 0.344852239587 | false |
| 100 | 9.89179256092 | 2.3798515067 | 2.96753776828 | 0.344852239587 | false |
| 10000 | 9.89161448599 | 2.37987968063 | 2.9674843458 | 0.34485215703 | false |
| 1000000 | 9.89337197345 | 2.37944423544 | 2.96801159204 | 0.344852555389 | false |
