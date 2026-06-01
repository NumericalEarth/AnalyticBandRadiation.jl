# Reduced ecCKD Retained Quadrature Linearized Optimizer

Status: **quadrature_linearized_optimizer_rejected**

| Field | Value |
|---|---:|
| Base mode | retained_topology |
| Residual mode | all_shortwave |
| Basis count | 15 |
| Probe step | 0.00390625 |
| Max logit delta | 0.01 |
| Base objective | 7.05452306522 |
| Base TOA forcing | 2.11635691956 W m^-2 |
| Base surface forcing | 2.0297810592 W m^-2 |
| Best ridge lambda | 10000 |
| Best exact objective | 8.92989013811 |
| Best objective reduction | -1.8753670729 |
| Best TOA forcing | 2.20185456657 W m^-2 |
| Best surface forcing | 2.67896704143 W m^-2 |
| Any Pareto-safe | false |
| Accepted | false |

This diagnostic fits a linearized all-logit shortwave quadrature-weight
update on the retained current base, then exact-evaluates each ridge
candidate and accepts only strict Pareto-safe TOA/surface improvements.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Pareto-safe |
|---:|---:|---:|---:|---:|
| 1e-08 | 10.1749754866 | 2.15572839315 | 3.05249264597 | false |
| 1e-06 | 10.1749754902 | 2.15572839322 | 3.05249264705 | false |
| 0.0001 | 10.1749758467 | 2.15572840083 | 3.05249275402 | false |
| 0.01 | 10.1750114961 | 2.15572916161 | 3.05250344883 | false |
| 1 | 10.1784837704 | 2.15580326204 | 3.05354513113 | false |
| 100 | 10.3299262001 | 2.15902135989 | 3.09897786003 | false |
| 10000 | 8.92989013811 | 2.20185456657 | 2.67896704143 | false |
| 1000000 | 9.33736791289 | 2.19266198126 | 2.80121037387 | false |
| 100000000 | 9.40879112697 | 2.215641526 | 2.82263733809 | false |
