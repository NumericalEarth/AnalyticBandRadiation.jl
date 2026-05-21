# Reduced ecCKD Retained Current Component-Scale Optimizer 3

Status: **current_component_scale_optimizer3_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer2_composed |
| Basis | third_pass_per_gpoint_static_h2o_rayleigh_component_scales |
| Basis count | 48 |
| Predecessor delta counts | 48, 48 |
| Probe step | 0.0009765625 |
| Max log scale | 0.00048828125 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.00025 |
| Base objective | 6.96730155579 |
| Base TOA forcing | 2.09019046674 W m^-2 |
| Base surface forcing | 2.01348577884 W m^-2 |
| Base heating RMSE | 0.346742121541 K day^-1 |
| Best objective | 6.96625146691 |
| Best objective reduction | 0.00105008888506 |
| Best TOA forcing | 2.08987544007 W m^-2 |
| Best surface forcing | 2.01423282621 W m^-2 |
| Best heating RMSE | 0.346419533777 K day^-1 |
| Accepted | true |
| Accepted ridge lambda | 10000 |
| Accepted objective | 6.96625146691 |
| Accepted TOA forcing | 2.08987544007 W m^-2 |
| Accepted surface forcing | 2.01423282621 W m^-2 |
| Accepted heating RMSE | 0.346419533777 K day^-1 |

This diagnostic composes the first two accepted component-scale moves,
then tests whether the same low-rank parameter family still has a
smaller cap-safe improving direction on the updated current base.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 6.96830564582 | 2.09049169375 | 2.01387953055 | 0.346331005641 | false |
| 1e-06 | 6.96709759444 | 2.09012927833 | 2.01444365268 | 0.346456063068 | false |
| 0.0001 | 6.96659598711 | 2.08997879613 | 2.01476168772 | 0.346446407012 | true |
| 0.01 | 6.96709444969 | 2.09012833491 | 2.01439656684 | 0.346457642211 | false |
| 1 | 6.96709036521 | 2.09012710956 | 2.01438691775 | 0.34645722025 | false |
| 100 | 6.96827155907 | 2.09048146772 | 2.01416197251 | 0.346474247201 | false |
| 10000 | 6.96625146691 | 2.08987544007 | 2.01423282621 | 0.346419533777 | true |
| 1000000 | 6.98337622519 | 2.09501286756 | 2.01520773454 | 0.34658974527 | false |
