# Reduced ecCKD Leave-One-Out Single Table-Move Scan

Status: **single_table_move_rejected**

| Field | Value |
|---|---:|
| Omitted SW g-point | 23 |
| Model | 32x31 |
| Candidate count | 8 |
| Exact rows | 48 |
| Pressure-index min | 1 |
| Pressure-index max | 8 |
| Boundary cap | 0.3 W m^-2 |
| Base objective | 1.62442314379 |
| Base boundary forcing | 0.0751392004466 W m^-2 |
| Base heating RMSE | 0.0812211571893 K day^-1 |
| Best objective | 12.4447761537 |
| Best objective reduction | -10.8203530099 |
| Best boundary forcing | 3.7334328461 W m^-2 |
| Best heating RMSE | 0.0932205723437 K day^-1 |
| Accepted | false |
| Accepted objective | 1.62442314379 |
| Accepted boundary forcing | 0.0751392004466 W m^-2 |
| Accepted heating RMSE | 0.0812211571893 K day^-1 |

This diagnostic exact-evaluates individual local active table-entry moves
on the objective-best 32x31 leave-one-out row by default. The default
pressure-index window targets the upper atmosphere where the g23
heating residual is localized. It is a cheap check for whether the
local table basis has any one-step descent direction before trying
another combined ridge step.

## Best Exact Rows

| Rank | Candidate | Log scale | Objective | Reduction | Boundary forcing | Heating RMSE | Accepted |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `component=static_absorption, local_gpoint_index=25, gpoint=26, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.001953125 | 12.4447761537 | -10.8203530099 | 3.7334328461 | 0.0932205723437 | false |
| 2 | `component=static_absorption, local_gpoint_index=1, gpoint=1, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.001953125 | 12.4447763314 | -10.8203531876 | 3.73343289941 | 0.093221049142 | false |
| 3 | `component=static_absorption, local_gpoint_index=25, gpoint=26, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.00048828125 | 12.4447763802 | -10.8203532364 | 3.73343291407 | 0.0932211656683 | false |
| 4 | `component=static_absorption, local_gpoint_index=19, gpoint=19, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.001953125 | 12.444776414 | -10.8203532702 | 3.7334329242 | 0.0932041892865 | false |
| 5 | `component=static_absorption, local_gpoint_index=1, gpoint=1, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.00048828125 | 12.4447764246 | -10.8203532808 | 3.73343292738 | 0.0932212847988 | false |
| 6 | `component=static_absorption, local_gpoint_index=2, gpoint=2, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.001953125 | 12.4447764271 | -10.8203532833 | 3.73343292812 | 0.0932203340375 | false |
| 7 | `component=static_absorption, local_gpoint_index=25, gpoint=26, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.0001220703125 | 12.4447764368 | -10.820353293 | 3.73343293104 | 0.0932213138857 | false |
| 8 | `component=static_absorption, local_gpoint_index=17, gpoint=17, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.001953125 | 12.4447764388 | -10.820353295 | 3.73343293164 | 0.0932206514272 | false |
| 9 | `component=static_absorption, local_gpoint_index=19, gpoint=19, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.00048828125 | 12.4447764452 | -10.8203533015 | 3.73343293357 | 0.0932170648039 | false |
| 10 | `component=static_absorption, local_gpoint_index=1, gpoint=1, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.0001220703125 | 12.4447764479 | -10.8203533041 | 3.73343293437 | 0.0932213436618 | false |
| 11 | `component=static_absorption, local_gpoint_index=2, gpoint=2, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.00048828125 | 12.4447764485 | -10.8203533047 | 3.73343293455 | 0.0932211061275 | false |
| 12 | `component=static_absorption, local_gpoint_index=29, gpoint=30, gas_index=1, pressure_index=3, temperature_index=4, h2o_index=0` | 0.001953125 | 12.4447764507 | -10.8203533069 | 3.73343293521 | 0.0932205445061 | false |
