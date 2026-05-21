# ecRad Clear Aerosol Optics Gap

This diagnostic compares generated clear gas+aerosol optical properties against ecRad's saved all-sky clear optical properties.

| Candidate kind | Variable | RMSE | Max abs | Mean bias | Mean abs | Reference mean | Candidate mean |
|---|---|---:|---:|---:|---:|---:|---:|
| `clear_gas_aerosol` | `od_lw` | 0.000876825556062 | 0.0128478864502 | 0.000228964146162 | 0.000278890961814 | 1.70518589148 | 1.70541485563 |
| `clear_gas_aerosol` | `od_sw` | 7.55486296587e-05 | 0.00133862378566 | 1.66520079124e-05 | 1.67099665583e-05 | 0.130029342962 | 0.13004599497 |
| `clear_gas_aerosol` | `ssa_sw` | 1.14757470859e-05 | 4.39661454027e-05 | 3.17786982342e-06 | 6.87575542339e-06 | 0.329059882431 | 0.329063060301 |
| `clear_gas_aerosol` | `asymmetry_sw` | 1.08435120764e-05 | 2.05960276292e-05 | -9.38797101172e-06 | 9.38797101172e-06 | 0.214922051934 | 0.214912663963 |

## Candidate Configuration

| Environment variable | Value |
|---|---|
| `RH_AEROSOL_DELTA_EDDINGTON_SCALE` | `true` |
| `RH_AEROSOL_OPTICS` | `true` |
| `RH_CANDIDATE_GAS_OPTICS` | `official_ecckd` |
| `RH_CLOUD_EFFECTIVE_RADIUS_OPTICS` | `false` |
| `RH_CLOUD_SCATTERING_TABLE_OPTICS` | `false` |
| `RH_IFS_AEROSOL_LAYER_RELATIVE_HUMIDITY` | `true` |
| `RH_IFS_AEROSOL_TABLE_OPTICS` | `true` |
