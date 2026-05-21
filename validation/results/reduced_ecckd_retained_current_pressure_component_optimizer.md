# Reduced ecCKD Retained Current Pressure-Component Optimizer

Status: **current_pressure_component_optimizer_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer4_composed |
| Basis | per_gpoint_pressure_band_static_h2o_component_scales |
| Basis count | 128 |
| Pressure bands | 4 |
| Probe step | 0.0009765625 |
| Max log scale | 0.001953125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.96570561645 |
| Base TOA forcing | 2.08971168493 W m^-2 |
| Base surface forcing | 2.01462006099 W m^-2 |
| Base heating RMSE | 0.346248776779 K day^-1 |
| Best objective | 6.9558996605 |
| Best objective reduction | 0.00980595595015 |
| Best TOA forcing | 2.08676989815 W m^-2 |
| Best surface forcing | 2.01383542392 W m^-2 |
| Best heating RMSE | 0.343997633773 K day^-1 |
| Accepted | true |
| Accepted ridge lambda | 1 |
| Accepted objective | 6.9558996605 |
| Accepted TOA forcing | 2.08676989815 W m^-2 |
| Accepted surface forcing | 2.01383542392 W m^-2 |
| Accepted heating RMSE | 0.343997633773 K day^-1 |
| Accepted moves | 119 |

This diagnostic composes the current four accepted component-scale
moves, then solves a pressure-band static/H2O component-scale basis
against the coupled heating-profile and boundary residual.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 7.00282502518 | 2.10084750756 | 2.00519074936 | 0.344785004704 | false |
| 1e-06 | 7.0028181366 | 2.10084544098 | 2.00519306427 | 0.345287180234 | false |
| 0.0001 | 7.00228572139 | 2.10068571642 | 2.00516155782 | 0.345080371864 | false |
| 0.01 | 7.00887280893 | 2.10266184268 | 2.00236695398 | 0.345455342844 | false |
| 1 | 6.9558996605 | 2.08676989815 | 2.01383542392 | 0.343997633773 | true |
| 100 | 6.98496490542 | 2.09548947163 | 2.0091952765 | 0.344462927639 | false |
| 10000 | 6.98285883898 | 2.09485765169 | 2.00755785349 | 0.343754819317 | false |
| 1000000 | 6.9649561678 | 2.08948685034 | 2.01431040285 | 0.344833293117 | false |
