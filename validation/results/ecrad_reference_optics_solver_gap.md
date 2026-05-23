# ecRad Reference-Optics Solver Gap

This diagnostic runs RadiativeHeating shortwave solvers using ecRad's saved all-sky shortwave optical properties. It isolates solver/source treatment from gas, cloud, and aerosol optical-property generation.

- Reference case: `validation/reference/ecrad/ecckd_all_sky_tropical_column.nc`
- ecRad properties: `validation/external/ecrad/test/ifs/radiative_properties_ecckd_tc.nc`
- ecRad output: `validation/external/ecrad/test/ifs/ecrad_meridian_ecckd_tc_out_REFERENCE.nc`

| Mode | SW up RMSE | SW down RMSE | Ref TOA net max abs | Ref surface net max abs | Output TOA net max abs | Output surface net max abs | Ref-output TOA max abs | Ref-output surface max abs | Clear direct max abs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `clear` | 4.00416960071e-06 | 3.83651624354e-05 | 7.00168743606e-06 | 2.64497741682e-05 | 5.27780546236e-05 | 4.80285355025e-05 | 5.34057617188e-05 | 3.0517578125e-05 | 6.09635194451e-05 |
| `grid_mean` | 100.850734965 | 56.2704550017 | 179.679290036 | 191.58353217 | 179.679259518 | 191.583537893 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `matrix_alpha` | 22.2642578063 | 15.0373640217 | 50.7779322971 | 48.5444545901 | 50.7779322971 | 48.544451729 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p0` | 13.0703922484 | 14.390447323 | 48.4759020663 | 44.1096923503 | 48.4759020663 | 44.1096875819 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p1` | 4.56883598055 | 4.49597394815 | 14.9511681022 | 13.3573949939 | 14.9511681022 | 13.3573902255 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p2` | 1.01811280215e-05 | 3.69130541194e-05 | 1.42029895187e-05 | 2.57961164607e-05 | 5.31249474989e-05 | 2.57961164607e-05 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p4` | 4.77708839394 | 4.15695565752 | 12.6545759312 | 11.1213983007 | 12.6545759312 | 11.121403069 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p8` | 8.91993034858 | 7.31940730999 | 22.8566633738 | 22.1126611668 | 22.8566633738 | 22.1126583057 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p12` | 10.820617469 | 8.61932942836 | 28.2181726774 | 27.2603471435 | 28.2181726774 | 27.2603442825 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
| `tripleclouds_alpha_p16` | 11.8951715532 | 9.30043744737 | 31.1488112397 | 30.0758785764 | 31.1488112397 | 30.0758757154 | 4.57763671875e-05 | 2.67028808594e-05 | n/a |
