# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer4_composed |
| Basis | per_gpoint_pressure_band_static_h2o_component_scales_plus_rayleigh |
| Include Rayleigh | true |
| Probe step | 0.0009765625 |
| Max log scale | 0.001953125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.96570561645 |
| Base TOA forcing | 2.08971168493 W m^-2 |
| Base surface forcing | 2.01462006099 W m^-2 |
| Base heating RMSE | 0.346248776779 K day^-1 |
| Selected partition | index |
| Selected pressure bands | 4 |
| Selected objective | 6.95607873366 |
| Selected TOA forcing | 2.0868236201 W m^-2 |
| Selected surface forcing | 2.02048400251 W m^-2 |
| Selected heating RMSE | 0.344138006611 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| index | true | 4 | 144 | 6.95607873366 | 2.0868236201 | 2.02048400251 | 0.344138006611 | true | 135 |
| index | true | 8 | 272 | 6.96570561645 | 2.08971168493 | 2.01462006099 | 0.346248776779 | false | 0 |
| log_pressure | true | 4 | 144 | 6.95732835204 | 2.08719850561 | 2.02009493142 | 0.343947778422 | true | 135 |
| log_pressure | true | 8 | 272 | 6.96570561645 | 2.08971168493 | 2.01462006099 | 0.346248776779 | false | 0 |
