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
| Max log scale | 0.0009765625 |
| Surface cap | 2.022 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.0005 |
| Base objective | 6.96570561645 |
| Base TOA forcing | 2.08971168493 W m^-2 |
| Base surface forcing | 2.01462006099 W m^-2 |
| Base heating RMSE | 0.346248776779 K day^-1 |
| Selected partition | log_pressure |
| Selected pressure bands | 4 |
| Selected objective | 6.9473445941 |
| Selected TOA forcing | 2.08420337823 W m^-2 |
| Selected surface forcing | 2.01957501177 W m^-2 |
| Selected heating RMSE | 0.34513970659 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Static gas split | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| log_pressure | false | true | 4 | 576 | 6.9473445941 | 2.08420337823 | 2.01957501177 | 0.34513970659 | true | 312 |
