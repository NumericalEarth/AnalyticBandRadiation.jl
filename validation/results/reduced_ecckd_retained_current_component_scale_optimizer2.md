# Reduced ecCKD Retained Current Component-Scale Optimizer 2

Status: **current_component_scale_optimizer2_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer_composed |
| Basis | second_pass_per_gpoint_static_h2o_rayleigh_component_scales |
| Basis count | 48 |
| First-pass delta count | 48 |
| Probe step | 0.0009765625 |
| Max log scale | 0.0009765625 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.0005 |
| Base objective | 6.96933451127 |
| Base TOA forcing | 2.09080035338 W m^-2 |
| Base surface forcing | 2.01498794456 W m^-2 |
| Base heating RMSE | 0.347361493009 K day^-1 |
| Best objective | 6.96730155579 |
| Best objective reduction | 0.00203295547702 |
| Best TOA forcing | 2.09019046674 W m^-2 |
| Best surface forcing | 2.01348577884 W m^-2 |
| Best heating RMSE | 0.346742121541 K day^-1 |
| Accepted | true |
| Accepted ridge lambda | 10000 |
| Accepted objective | 6.96730155579 |
| Accepted TOA forcing | 2.09019046674 W m^-2 |
| Accepted surface forcing | 2.01348577884 W m^-2 |
| Accepted heating RMSE | 0.346742121541 K day^-1 |

This diagnostic composes the accepted first component-scale move, then
tests whether the same low-rank parameter family has another cap-safe
improving direction on the updated current base.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 6.97134485039 | 2.09140345512 | 2.02563069571 | 0.346542556489 | false |
| 1e-06 | 6.96892588594 | 2.09067776578 | 2.02358100339 | 0.346791836809 | false |
| 0.0001 | 6.96792102511 | 2.09037630753 | 2.0235990981 | 0.346772451765 | true |
| 0.01 | 6.968919462 | 2.0906758386 | 2.02422833773 | 0.346794971663 | false |
| 1 | 6.96891134348 | 2.09067340304 | 2.02467039007 | 0.346794139816 | false |
| 100 | 6.9712791015 | 2.09138373045 | 2.01336843656 | 0.346824274951 | false |
| 10000 | 6.96730155579 | 2.09019046674 | 2.01348577884 | 0.346742121541 | true |
| 1000000 | 6.99807594198 | 2.0994227826 | 2.01666334192 | 0.347058637096 | false |
