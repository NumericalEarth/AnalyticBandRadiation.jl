# ecRad All-Sky Optics Gap

This diagnostic compares candidate optical properties against ecRad's saved all-sky radiative properties before radiative-transfer solving.

- Reference case: `validation/reference/ecrad/ecckd_32x64_all_sky_tropical_column.nc`
- ecRad properties: `validation/external/ecrad/test/ifs/radiative_properties.nc`

| Candidate kind | Variable | RMSE | Max abs | Mean bias | Mean abs | Reference mean | Candidate mean |
|---|---|---:|---:|---:|---:|---:|---:|
| `clear_region_current` | `od_lw` | 0.000905458790644 | 0.0129436451471 | 0.000275037795067 | 0.000277953983832 | 1.70518589148 | 1.70546092928 |
| `clear_region_current` | `od_sw` | 0.000169568051809 | 0.00211905626518 | 2.65310639835e-05 | 3.38168164462e-05 | 0.580674604003 | 0.580701135067 |
| `clear_region_current` | `ssa_sw` | 0.00196772932312 | 0.0510614898224 | -0.000463154552284 | 0.000509061761108 | 0.29246840701 | 0.292005252458 |
| `clear_region_current` | `asymmetry_sw` | 0.0100876046307 | 0.0532396749776 | -0.0040275168305 | 0.00426259761099 | 0.205867224332 | 0.201839707502 |
| `cloudy_region_delta_scaled` | `od_lw_cloud` | 1.22097248677 | 16.1661612906 | 0.283800617846 | 0.283800617846 | 0.512200058516 | 0.796000676362 |
| `cloudy_region_delta_scaled` | `ssa_lw_cloud` | 0.149091955548 | 0.490818006552 | 0.0636352559752 | 0.0636352559752 | 0.0410393158127 | 0.104674571788 |
| `cloudy_region_delta_scaled` | `asymmetry_lw_cloud` | 0.180972508112 | 0.478504217909 | 0.0807217288524 | 0.0807217288524 | 0.0944003760312 | 0.175122104884 |
| `cloudy_region_delta_scaled` | `od_sw_cloud` | 0.0657967006849 | 1.75544145773 | 0.00691016154387 | 0.00738528361851 | 0.258591710196 | 0.26550187174 |
| `cloudy_region_delta_scaled` | `ssa_sw_cloud` | 0.00647557156045 | 0.0902435067871 | -0.000121642650045 | 0.00128960161751 | 0.164386830278 | 0.164265187628 |
| `cloudy_region_delta_scaled` | `asymmetry_sw_cloud` | 0.00156598537057 | 0.0270432701402 | -0.000182887724493 | 0.000237818081822 | 0.0952936659519 | 0.0951107782274 |
| `cloudy_region_delta_scaled` | `od_lw+od_lw_cloud` | 1.22109993023 | 16.1661612254 | 0.284075655641 | 0.284077394595 | 2.21738595 | 2.50146160564 |
| `cloudy_region_delta_scaled` | `od_sw+od_sw_cloud` | 0.0658072946626 | 1.75554005586 | 0.00693669260786 | 0.00741611990026 | 0.839266314198 | 0.846203006806 |
| `cloudy_region_delta_scaled` | `combined_ssa_sw` | 0.00561894802282 | 0.0900858990526 | -0.000269353981514 | 0.00129599436699 | 0.36998143424 | 0.369712080259 |
| `cloudy_region_delta_scaled` | `combined_asymmetry_sw` | 0.00997744526199 | 0.0532396749776 | -0.00385103311996 | 0.00407133443283 | 0.243411927418 | 0.239560894298 |
| `boundary_materialized` | `incoming_sw` | 0 | 0 | 0 | 0 | 22.0044536591 | 22.0044536591 |
| `boundary_materialized` | `sw_albedo` | 0 | 0 | 0 | 0 | 0.0978319282844 | 0.0978319282844 |
| `boundary_materialized` | `sw_albedo_direct` | 0 | 0 | 0 | 0 | 0.0809983684143 | 0.0809983684143 |

## Candidate Configuration

| Environment variable | Value |
|---|---|
| `RH_AEROSOL_OPTICS` | `true` |
| `RH_CANDIDATE_GAS_OPTICS` | `official_ecckd` |
| `RH_CLOUD_EFFECTIVE_RADIUS_OPTICS` | `true` |
| `RH_CLOUD_FRACTION_EXPONENT` | `1.0` |
| `RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION` | `false` |
| `RH_CLOUD_OVERLAP_RULE` | `tripleclouds_alpha` |
| `RH_CLOUD_OVERLAP_SHORTWAVE` | `true` |
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
| `RH_IFS_AEROSOL_TABLE_OPTICS` | `true` |
| `RH_LIQUID_CLOUD_EFFECTIVE_RADIUS_SW_SCALE` | `1.0` |
| `RH_LIQUID_CLOUD_LW_MASS_ABSORPTION` | `100` |
| `RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY` | `0.85` |
| `RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO` | `1.0` |
| `RH_ECCKD_SW_PATH` | `validation/external/ecrad/data/ecckd-1.2_sw_climate_window-64b_ckd-definition.nc` |
| `RH_CLOUD_SW_MAPPING_PATH` | `validation/external/ecrad/data/ecckd-1.2_sw_climate_window-64b_ckd-definition.nc` |
| `RH_LW_CLOUD_SCATTERING` | `true` |
| `RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS` | `true` |
| `RH_CLOUD_INHOM_OVERLAP_EXPONENT` | `2.0` |
| `RH_CLOUD_OVERLAP_LONGWAVE` | `true` |
| `RH_CLOUD_OVERLAP_LONGWAVE_RULE` | `tripleclouds_alpha` |
| `RH_ALL_SKY_OPTICS_REFERENCE` | `validation/reference/ecrad/ecckd_32x64_all_sky_tropical_column.nc` |
| `RH_ALL_SKY_OPTICS_PROPERTIES` | `validation/external/ecrad/test/ifs/radiative_properties.nc` |
| `RH_ALL_SKY_OPTICS_OUTPUT_BASENAME` | `ecrad_all_sky_optics_gap_32x64_gate` |
