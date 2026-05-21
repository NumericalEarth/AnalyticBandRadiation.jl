# ecRad Reference-Optics Solver Gap

This diagnostic runs RadiativeHeating shortwave solvers using ecRad's saved all-sky shortwave optical properties. It isolates solver/source treatment from gas, cloud, and aerosol optical-property generation.

| Mode | SW up RMSE | SW down RMSE | TOA net max abs | Surface net max abs | TOA net bias | Surface net bias | Clear direct max abs |
|---|---:|---:|---:|---:|---:|---:|---:|
| `clear` | 4.00416960071e-06 | 3.83651624354e-05 | 7.00168743606e-06 | 2.64497741682e-05 | -6.17966509253e-07 | 1.09977016791e-06 | 6.09635194451e-05 |
| `grid_mean` | 100.850734965 | 56.2704550017 | 179.679290036 | 191.58353217 | -79.941222525 | -84.9522253174 | n/a |
| `matrix_alpha` | 22.2642578063 | 15.0373640217 | 50.7779322971 | 48.5444545901 | -18.8064335377 | -18.9913659036 | n/a |
| `tripleclouds_alpha_p0` | 13.0703922484 | 14.390447323 | 48.4759020663 | 44.1096923503 | 9.2052411121 | 8.95941727081 | n/a |
| `tripleclouds_alpha_p1` | 4.56883598055 | 4.49597394815 | 14.9511681022 | 13.3573949939 | 3.26666596151 | 3.16618147191 | n/a |
| `tripleclouds_alpha_p2` | 1.01811280215e-05 | 3.69130541194e-05 | 1.42029895187e-05 | 2.57961164607e-05 | -5.77021522759e-06 | 3.21874293547e-06 | n/a |
| `tripleclouds_alpha_p4` | 4.77708839394 | 4.15695565752 | 12.6545759312 | 11.1213983007 | -3.53254076656 | -3.4330284528 | n/a |
| `tripleclouds_alpha_p8` | 8.91993034858 | 7.31940730999 | 22.8566633738 | 22.1126611668 | -6.75437981256 | -6.59404447985 | n/a |
| `tripleclouds_alpha_p12` | 10.820617469 | 8.61932942836 | 28.2181726774 | 27.2603471435 | -8.32858072033 | -8.1605112413 | n/a |
| `tripleclouds_alpha_p16` | 11.8951715532 | 9.30043744737 | 31.1488112397 | 30.0758785764 | -9.26371501064 | -9.1024582155 | n/a |
