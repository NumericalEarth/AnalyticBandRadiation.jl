# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer4_composed |
| Basis | per_gpoint_pressure_band_static_gas_h2o_component_scales |
| Include Rayleigh | false |
| Static gas split | true |
| Probe step | 0.0009765625 |
| Max log scale | 0.001953125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.96570561645 |
| Base TOA forcing | 2.08971168493 W m^-2 |
| Base surface forcing | 2.01462006099 W m^-2 |
| Base heating RMSE | 0.346248776779 K day^-1 |
| Selected partition | log_pressure |
| Selected pressure bands | 2 |
| Selected objective | 6.93636845447 |
| Selected TOA forcing | 2.08091053634 W m^-2 |
| Selected surface forcing | 2.02313898324 W m^-2 |
| Selected heating RMSE | 0.344896146818 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Static gas split | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| log_pressure | false | true | 2 | 288 | 6.93636845447 | 2.08091053634 | 2.02313898324 | 0.344896146818 | true | 158 |
