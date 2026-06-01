# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer4_composed |
| Basis | per_gpoint_pressure_band_static_h2o_component_scales |
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
| Selected pressure bands | 8 |
| Selected objective | 6.93440630944 |
| Selected TOA forcing | 2.08032189283 W m^-2 |
| Selected surface forcing | 2.02314211275 W m^-2 |
| Selected heating RMSE | 0.345536576619 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| index | 2 | 64 | 6.9512588103 | 2.08537764309 | 2.01571378725 | 0.345036307871 | true | 60 |
| index | 4 | 128 | 6.9558996605 | 2.08676989815 | 2.01383542392 | 0.343997633773 | true | 119 |
| index | 8 | 256 | 6.9365355379 | 2.08096066137 | 2.02248785744 | 0.344227964572 | true | 238 |
| index | 16 | 512 | 6.95995780167 | 2.0879873405 | 2.01505360681 | 0.344492267568 | true | 474 |
| log_pressure | 2 | 64 | 6.95043162612 | 2.08512948784 | 2.01596431868 | 0.344055128559 | true | 60 |
| log_pressure | 4 | 128 | 6.95430043407 | 2.08629013022 | 2.01452212417 | 0.343995548869 | true | 119 |
| log_pressure | 8 | 256 | 6.93440630944 | 2.08032189283 | 2.02314211275 | 0.345536576619 | true | 238 |
| log_pressure | 16 | 512 | 6.94611808289 | 2.08383542487 | 2.02085868528 | 0.344817750125 | true | 474 |
