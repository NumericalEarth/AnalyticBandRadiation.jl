# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_gas_pressure_component_scan_composed |
| Basis | per_gpoint_pressure_band_static_gas_h2o_component_scales |
| Include Rayleigh | false |
| Static gas split | true |
| Heating residual weight | 1 |
| Boundary residual weight | 1 |
| Probe step | 0.0009765625 |
| Max log scale | 0.00048828125 |
| Surface cap | 2.028 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.00025 |
| Base objective | 6.92896358191 |
| Base TOA forcing | 2.07868907457 W m^-2 |
| Base surface forcing | 2.02453550196 W m^-2 |
| Base heating RMSE | 0.344039253986 K day^-1 |
| Selected partition | log_pressure |
| Selected pressure bands | 4 |
| Selected objective | 6.919764946 |
| Selected TOA forcing | 2.0759294838 W m^-2 |
| Selected surface forcing | 2.02701815786 W m^-2 |
| Selected heating RMSE | 0.343492243875 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Static gas split | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| log_pressure | false | true | 4 | 576 | 6.919764946 | 2.0759294838 | 2.02701815786 | 0.343492243875 | true | 312 |
