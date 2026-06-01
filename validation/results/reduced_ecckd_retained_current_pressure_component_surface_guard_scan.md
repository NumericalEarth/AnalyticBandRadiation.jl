# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer4_composed |
| Basis | per_gpoint_pressure_band_static_h2o_component_scales |
| Include Rayleigh | false |
| Probe step | 0.0009765625 |
| Max log scale | 0.001953125 |
| Surface cap | 2.02 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.001 |
| Base objective | 6.96570561645 |
| Base TOA forcing | 2.08971168493 W m^-2 |
| Base surface forcing | 2.01462006099 W m^-2 |
| Base heating RMSE | 0.346248776779 K day^-1 |
| Selected partition | log_pressure |
| Selected pressure bands | 4 |
| Selected objective | 6.95430043407 |
| Selected TOA forcing | 2.08629013022 W m^-2 |
| Selected surface forcing | 2.01452212417 W m^-2 |
| Selected heating RMSE | 0.343995548869 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| index | false | 4 | 128 | 6.9558996605 | 2.08676989815 | 2.01383542392 | 0.343997633773 | true | 119 |
| index | false | 8 | 256 | 6.96570561645 | 2.08971168493 | 2.01462006099 | 0.346248776779 | false | 0 |
| log_pressure | false | 4 | 128 | 6.95430043407 | 2.08629013022 | 2.01452212417 | 0.343995548869 | true | 119 |
| log_pressure | false | 8 | 256 | 6.96570561645 | 2.08971168493 | 2.01462006099 | 0.346248776779 | false | 0 |
