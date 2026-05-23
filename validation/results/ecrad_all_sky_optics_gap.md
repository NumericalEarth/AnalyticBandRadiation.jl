# ecRad All-Sky Optics Gap

This diagnostic compares candidate optical properties against ecRad's saved all-sky radiative properties before radiative-transfer solving.

- Reference case: `validation/reference/ecrad/ecckd_all_sky_tropical_column.nc`
- ecRad properties: `validation/external/ecrad/test/ifs/radiative_properties_ecckd_tc.nc`

| Candidate kind | Variable | RMSE | Max abs | Mean bias | Mean abs | Reference mean | Candidate mean |
|---|---|---:|---:|---:|---:|---:|---:|
| `clear_region_current` | `od_lw` | 0.000879836629385 | 0.0123799379257 | 0.000100906121824 | 0.000333790186911 | 1.70518589148 | 1.7052867976 |
| `clear_region_current` | `od_sw` | 0.000993550204802 | 0.00590667874316 | -0.000499946495323 | 0.000515142452707 | 0.130029342962 | 0.129529396467 |
| `clear_region_current` | `ssa_sw` | 0.121462193171 | 0.785962885946 | -0.051650669114 | 0.0594167640114 | 0.329059882431 | 0.277409213317 |
| `clear_region_current` | `asymmetry_sw` | 0.253536812807 | 0.464174984389 | -0.214922051934 | 0.214922051934 | 0.214922051934 | 0 |
| `cloudy_region_delta_scaled` | `od_lw_cloud` | 1.22097248677 | 16.1661612906 | 0.283800617846 | 0.283800617846 | 0.512200058516 | 0.796000676362 |
| `cloudy_region_delta_scaled` | `ssa_lw_cloud` | 0.149091955548 | 0.490818006552 | 0.0636352559752 | 0.0636352559752 | 0.0410393158127 | 0.104674571788 |
| `cloudy_region_delta_scaled` | `asymmetry_lw_cloud` | 0.180972508112 | 0.478504217909 | 0.0807217288524 | 0.0807217288524 | 0.0944003760312 | 0.175122104884 |
| `cloudy_region_delta_scaled` | `od_sw_cloud` | 6.28953105577e-05 | 0.00136237868001 | 7.65853442876e-06 | 9.32993948313e-06 | 0.255843425011 | 0.255851083545 |
| `cloudy_region_delta_scaled` | `ssa_sw_cloud` | 7.71088598755e-06 | 0.000214932524265 | -1.40981514465e-06 | 1.43137537163e-06 | 0.173488713372 | 0.173487303557 |
| `cloudy_region_delta_scaled` | `asymmetry_sw_cloud` | 5.05537374548e-07 | 3.42462552154e-05 | 4.80197062469e-08 | 6.57417776473e-08 | 0.0950404614059 | 0.0950405094256 |
| `cloudy_region_delta_scaled` | `od_lw+od_lw_cloud` | 1.22102964154 | 16.1661128394 | 0.283901523968 | 0.284047330604 | 2.21738595 | 2.50128747396 |
| `cloudy_region_delta_scaled` | `od_sw+od_sw_cloud` | 0.00099128729058 | 0.00590668271613 | -0.000492287960894 | 0.000513036145906 | 0.385872767973 | 0.385380480012 |
| `cloudy_region_delta_scaled` | `combined_ssa_sw` | 0.106099653537 | 0.785962885946 | -0.0402074426108 | 0.0462873204548 | 0.416353085541 | 0.37614564293 |
| `cloudy_region_delta_scaled` | `combined_asymmetry_sw` | 0.211834879183 | 0.464174984389 | -0.157602882686 | 0.157692706575 | 0.250728665661 | 0.0931257829755 |
| `boundary_materialized` | `incoming_sw` | 17.1081617253 | 82.0313539575 | -9.21656971893 | 9.21656971893 | 44.0089073181 | 34.7923375992 |
| `boundary_materialized` | `sw_albedo` | 0 | 0 | 0 | 0 | 0.100239364205 | 0.100239364205 |
| `boundary_materialized` | `sw_albedo_direct` | 0 | 0 | 0 | 0 | 0.0832147123029 | 0.0832147123029 |

## Candidate Configuration

| Environment variable | Value |
|---|---|
| `RH_AEROSOL_OPTICS` | `false` |
| `RH_CANDIDATE_GAS_OPTICS` | `official_ecckd` |
| `RH_CLOUD_EFFECTIVE_RADIUS_OPTICS` | `true` |
| `RH_CLOUD_FRACTION_EXPONENT` | `0.5` |
| `RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION` | `false` |
| `RH_CLOUD_OVERLAP_RULE` | `maximum` |
| `RH_CLOUD_OVERLAP_SHORTWAVE` | `false` |
| `RH_CLOUD_SCATTERING_DELTA_EDDINGTON_AVERAGE` | `true` |
| `RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE` | `true` |
| `RH_CLOUD_SCATTERING_MAPPING_METHOD` | `ecrad` |
| `RH_CLOUD_SCATTERING_SW_SSA_SCALE` | `1.0` |
| `RH_CLOUD_SCATTERING_TABLE_OPTICS` | `true` |
| `RH_CLOUD_SCATTERING_THICK_AVERAGING` | `true` |
| `RH_ICE_CLOUD_EFFECTIVE_RADIUS_SW_SCALE` | `1.0` |
| `RH_ICE_CLOUD_LW_MASS_ABSORPTION` | `50` |
| `RH_ICE_CLOUD_SW_SCATTERING_ASYMMETRY` | `0.75` |
| `RH_ICE_CLOUD_SW_SINGLE_SCATTERING_ALBEDO` | `1.0` |
| `RH_IFS_AEROSOL_TABLE_OPTICS` | `false` |
| `RH_LIQUID_CLOUD_EFFECTIVE_RADIUS_SW_SCALE` | `1.0` |
| `RH_LIQUID_CLOUD_LW_MASS_ABSORPTION` | `100` |
| `RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY` | `0.85` |
| `RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO` | `1.0` |
