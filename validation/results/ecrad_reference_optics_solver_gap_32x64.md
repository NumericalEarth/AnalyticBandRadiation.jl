# ecRad Reference-Optics Solver Gap

This diagnostic runs RadiativeHeating shortwave solvers using ecRad's saved all-sky shortwave optical properties. It isolates solver/source treatment from gas, cloud, and aerosol optical-property generation.

- Reference case: `validation/reference/ecrad/ecckd_32x64_all_sky_tropical_column.nc`
- ecRad properties: `validation/external/ecrad/test/ifs/radiative_properties.nc`
- ecRad output: `validation/external/ecrad/test/ifs/ecrad_meridian_ecckd_32x64_all_sky_props_out.nc`

| Mode | SW up RMSE | SW down RMSE | Ref TOA net max abs | Ref surface net max abs | Output TOA net max abs | Output surface net max abs | Ref-output TOA max abs | Ref-output surface max abs | Clear direct max abs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `clear` | 4.17796906995e-06 | 3.81406437279e-05 | 5.93723370912e-06 | 3.91987397279e-05 | 3.81363654469e-05 | 4.26767421686e-05 | 3.81469726562e-05 | 3.0517578125e-05 | 6.10174265603e-05 |
| `grid_mean` | 100.862606143 | 56.0495373518 | 179.762920359 | 190.872992416 | 179.762920359 | 190.872986693 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `matrix_alpha` | 22.2741074785 | 15.0011605068 | 50.8938241721 | 48.5933699357 | 50.8938241721 | 48.593371843 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p0` | 12.9532014014 | 14.3703816778 | 47.998668362 | 44.048122106 | 47.998668362 | 44.0481182913 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p1` | 4.54180106923 | 4.49352060471 | 14.8402291529 | 13.3554732217 | 14.8402291529 | 13.3554694071 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p2` | 1.04877004832e-05 | 3.82551446446e-05 | 2.08785033919e-05 | 3.91987397279e-05 | 5.09760238856e-05 | 6.04556361168e-05 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p4` | 4.76342643824 | 4.15821994051 | 12.5875967386 | 11.1346554722 | 12.5875967386 | 11.1346592869 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p8` | 8.90623348943 | 7.32147317299 | 22.9095814619 | 22.1613553888 | 22.9095814619 | 22.1613572961 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p12` | 10.809983154 | 8.6204489493 | 28.2879703018 | 27.3144392088 | 28.2879703018 | 27.3144411162 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
| `tripleclouds_alpha_p16` | 11.8869028481 | 9.30038727276 | 31.2285281698 | 30.1317289689 | 31.2285281698 | 30.1317308763 | 4.57763671875e-05 | 3.0517578125e-05 | n/a |
