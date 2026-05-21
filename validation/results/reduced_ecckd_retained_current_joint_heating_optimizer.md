# Reduced ecCKD Retained Current Joint Heating Optimizer

Status: **current_joint_heating_optimizer_rejected**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Basis | quadrature_logits_plus_active_table_entries |
| Logit basis count | 15 |
| Table candidate count | 8 |
| Basis count | 23 |
| Logit probe step | 0.0009765625 |
| Table probe step | 0.00048828125 |
| Max logit delta | 0.00390625 |
| Max table log scale | 0.001953125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.97212754254 |
| Base TOA forcing | 2.09163826276 W m^-2 |
| Base surface forcing | 2.00694307171 W m^-2 |
| Base heating RMSE | 0.348575997036 K day^-1 |
| Best objective | 8.51798339957 |
| Best objective reduction | -1.54585585702 |
| Best TOA forcing | 2.49840725461 W m^-2 |
| Best surface forcing | 2.55539501987 W m^-2 |
| Best heating RMSE | 0.340399423623 K day^-1 |
| Accepted | false |
| Accepted ridge lambda | n/a |
| Accepted objective | 6.97212754254 |
| Accepted TOA forcing | 2.09163826276 W m^-2 |
| Accepted surface forcing | 2.00694307171 W m^-2 |
| Accepted heating RMSE | 0.348575997036 K day^-1 |
| Accepted table moves | 0 |

This diagnostic fits one local linear system with both shortwave
quadrature-logit columns and active table-entry columns. It tests
whether weight redistribution can offset the boundary regression seen
in table-only heating-profile residual fits while preserving the same
hard objective, TOA, surface, and heating-RMSE acceptance guards.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 8.97381544255 | 2.46778731542 | 2.69214463277 | 0.342036920497 | false |
| 1e-06 | 9.73221702526 | 2.38810506777 | 2.91966510758 | 0.339945174944 | false |
| 0.0001 | 8.88445526202 | 2.47061985932 | 2.6653365786 | 0.340331995186 | false |
| 0.01 | 8.8460249958 | 2.47354513524 | 2.65380749874 | 0.340339432064 | false |
| 1 | 8.51798339957 | 2.49840725461 | 2.55539501987 | 0.340399423623 | false |
| 100 | 9.73168820659 | 2.38817599022 | 2.91950646198 | 0.339944904901 | false |
| 10000 | 10.2936475441 | 2.23431676113 | 3.08809426323 | 0.343425469111 | false |
| 1000000 | 8.98637873849 | 2.34215302937 | 2.69591362155 | 0.343751369539 | false |
