# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_rejected**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_gas_pressure_component_continuation_scan_composed |
| Basis | per_gpoint_pressure_band_static_gas_h2o_component_scales |
| Include Rayleigh | false |
| Static gas split | true |
| Heating residual weight | 2 |
| Boundary residual weight | 1 |
| Probe step | 0.0009765625 |
| Max log scale | 0.000244140625 |
| Surface cap | 2.02702 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.000125 |
| Base objective | 6.919764946 |
| Base TOA forcing | 2.0759294838 W m^-2 |
| Base surface forcing | 2.02701815786 W m^-2 |
| Base heating RMSE | 0.343492243875 K day^-1 |
| Selected partition |  |
| Selected pressure bands | 0 |
| Selected objective | 6.919764946 |
| Selected TOA forcing | 2.0759294838 W m^-2 |
| Selected surface forcing | 2.02701815786 W m^-2 |
| Selected heating RMSE | 0.343492243875 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Static gas split | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| log_pressure | false | true | 4 | 576 | 6.919764946 | 2.0759294838 | 2.02701815786 | 0.343492243875 | false | 0 |
