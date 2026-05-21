# Reduced ecCKD Retained Current Component-Scale Optimizer 4

Status: **current_component_scale_optimizer4_improved**

| Field | Value |
|---|---:|
| Residual mode | heating_profile_boundary |
| Base mode | current_component_scale_optimizer3_composed |
| Basis | fourth_pass_per_gpoint_static_h2o_rayleigh_component_scales |
| Basis count | 48 |
| Predecessor delta counts | 48, 48, 48 |
| Probe step | 0.0009765625 |
| Max log scale | 0.000244140625 |
| Surface cap | 2.03 W m^-2 |
| TOA tolerance | 0 W m^-2 |
| Minimum objective reduction | 0.000125 |
| Base objective | 6.96625146691 |
| Base TOA forcing | 2.08987544007 W m^-2 |
| Base surface forcing | 2.01423282621 W m^-2 |
| Base heating RMSE | 0.346419533777 K day^-1 |
| Best objective | 6.96570561645 |
| Best objective reduction | 0.000545850463851 |
| Best TOA forcing | 2.08971168493 W m^-2 |
| Best surface forcing | 2.01462006099 W m^-2 |
| Best heating RMSE | 0.346248776779 K day^-1 |
| Accepted | true |
| Accepted ridge lambda | 10000 |
| Accepted objective | 6.96570561645 |
| Accepted TOA forcing | 2.08971168493 W m^-2 |
| Accepted surface forcing | 2.01462006099 W m^-2 |
| Accepted heating RMSE | 0.346248776779 K day^-1 |

This diagnostic composes the first three accepted component-scale moves,
then tests whether the same low-rank parameter family still has a
smaller cap-safe improving direction on the updated current base.

## Ridge Summary

| Ridge lambda | Objective | TOA forcing | Surface forcing | Heating RMSE | Accepted |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 6.96675322237 | 2.09002596671 | 2.01443026896 | 0.346213670148 | false |
| 1e-06 | 6.96614956132 | 2.0898448684 | 2.01471219275 | 0.346276246801 | false |
| 0.0001 | 6.96589896387 | 2.08976968916 | 2.01487110328 | 0.346271427972 | true |
| 0.01 | 6.96614800684 | 2.08984440205 | 2.01468863391 | 0.346277038579 | false |
| 1 | 6.96614595845 | 2.08984378753 | 2.01468381204 | 0.346276825901 | false |
| 100 | 6.96673641479 | 2.09002092444 | 2.0145713885 | 0.34628578969 | false |
| 10000 | 6.96570561645 | 2.08971168493 | 2.01462006099 | 0.346248776779 | true |
| 1000000 | 6.97521760499 | 2.0925652815 | 2.0148345329 | 0.346346590432 | false |
