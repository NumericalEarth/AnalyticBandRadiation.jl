# Reduced ecCKD Retained Current Pressure-Component Scan

Status: **current_pressure_component_scan_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_gas_pressure_component_continuation2_scan_composed |
| Basis | per_gpoint_pressure_band_static_gas_h2o_component_scales |
| Include Rayleigh | false |
| Static gas split | true |
| Heating residual weight | 8 |
| Boundary residual weight | 1 |
| Probe step | 0.0009765625 |
| Max log scale | 0.0001220703125 |
| Surface cap | 2.0270142 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 6.25e-05 |
| Base objective | 6.91802050742 |
| Base TOA forcing | 2.07540615223 W m^-2 |
| Base surface forcing | 2.02701415411 W m^-2 |
| Base heating RMSE | 0.343237488292 K day^-1 |
| Selected partition | log_pressure |
| Selected pressure bands | 4 |
| Selected objective | 6.91747223736 |
| Selected TOA forcing | 2.07524167121 W m^-2 |
| Selected surface forcing | 2.02695668832 W m^-2 |
| Selected heating RMSE | 0.343116151859 K day^-1 |

This diagnostic compares pressure-component bases before promoting any
new pressure-band partition into the canonical reduced row.

## Variants

| Partition | Rayleigh | Static gas split | Pressure bands | Basis count | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted | Moves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| log_pressure | false | true | 4 | 576 | 6.91747223736 | 2.07524167121 | 2.02695668832 | 0.343116151859 | true | 313 |
