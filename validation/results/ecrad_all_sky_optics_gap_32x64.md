# ecRad All-Sky Optics Gap

This diagnostic compares candidate optical properties against ecRad's saved all-sky radiative properties before radiative-transfer solving.

- Reference case: `validation/reference/ecrad/ecckd_32x64_all_sky_tropical_column.nc`
- ecRad properties: `validation/external/ecrad/test/ifs/radiative_properties.nc`

| Candidate kind | Variable | RMSE | Max abs | Mean bias | Mean abs | Reference mean | Candidate mean |
|---|---|---:|---:|---:|---:|---:|---:|
| `clear_region_current` | `od_lw` | 0.000879836629385 | 0.0123799379257 | 0.000100906121824 | 0.000333790186911 | 1.70518589148 | 1.7052867976 |
| `clear_region_current` | `od_sw` | 0.0010618479152 | 0.00598724640426 | -0.000520998522495 | 0.000556629177587 | 0.580674604003 | 0.58015360548 |
| `clear_region_current` | `ssa_sw` | 0.099964356407 | 0.807057874438 | -0.0367335538063 | 0.0442113622658 | 0.29246840701 | 0.255734853204 |
| `clear_region_current` | `asymmetry_sw` | 0.25001091365 | 0.472405069614 | -0.205867224332 | 0.205867224332 | 0.205867224332 | 0 |
| `cloudy_region_delta_scaled` | `od_lw_cloud` | 1.22097248677 | 16.1661612906 | 0.283800617846 | 0.283800617846 | 0.512200058516 | 0.796000676362 |
| `cloudy_region_delta_scaled` | `ssa_lw_cloud` | 0.149091955548 | 0.490818006552 | 0.0636352559752 | 0.0636352559752 | 0.0410393158127 | 0.104674571788 |
| `cloudy_region_delta_scaled` | `asymmetry_lw_cloud` | 0.180972508112 | 0.478504217909 | 0.0807217288524 | 0.0807217288524 | 0.0944003760312 | 0.175122104884 |
| `cloudy_region_delta_scaled` | `od_sw_cloud` | 0.0657967006849 | 1.75544145773 | 0.00691016154387 | 0.00738528361851 | 0.258591710196 | 0.26550187174 |
| `cloudy_region_delta_scaled` | `ssa_sw_cloud` | 0.00647557156045 | 0.0902435067871 | -0.000121642650045 | 0.00128960161751 | 0.164386830278 | 0.164265187628 |
| `cloudy_region_delta_scaled` | `asymmetry_sw_cloud` | 0.00156598537057 | 0.0270432701402 | -0.000182887724493 | 0.000237818081822 | 0.0952936659519 | 0.0951107782274 |
| `cloudy_region_delta_scaled` | `od_lw+od_lw_cloud` | 1.22102964154 | 16.1661128394 | 0.283901523968 | 0.284047330604 | 2.21738595 | 2.50128747396 |
| `cloudy_region_delta_scaled` | `od_sw+od_sw_cloud` | 0.0657565286686 | 1.75543833319 | 0.00638916302138 | 0.00782493992001 | 0.839266314198 | 0.84565547722 |
| `cloudy_region_delta_scaled` | `combined_ssa_sw` | 0.0880324009694 | 0.807057874438 | -0.0289229334427 | 0.0359415016432 | 0.36998143424 | 0.341058500798 |
| `cloudy_region_delta_scaled` | `combined_asymmetry_sw` | 0.208603307371 | 0.472405069614 | -0.150830123403 | 0.150963554079 | 0.243411927418 | 0.0925818040154 |
| `boundary_materialized` | `incoming_sw` | 11.5323595963 | 82.0313546779 | -4.60828482966 | 4.60828482966 | 22.0044536591 | 17.3961688294 |
| `boundary_materialized` | `sw_albedo` | 0.0440527261062 | 0.288952410221 | -0.0219137080855 | 0.0219137080855 | 0.0978319282844 | 0.0759182201989 |
| `boundary_materialized` | `sw_albedo_direct` | 0.0380128835272 | 0.273333489895 | -0.00508014821532 | 0.0133205662147 | 0.0809983684143 | 0.0759182201989 |

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
| `RH_ECCKD_SW_PATH` | `validation/external/ecrad/data/ecckd-1.2_sw_climate_window-64b_ckd-definition.nc` |
| `RH_CLOUD_SW_MAPPING_PATH` | `validation/external/ecrad/data/ecckd-1.2_sw_climate_window-64b_ckd-definition.nc` |
| `RH_LW_CLOUD_SCATTERING` | `false` |
| `RH_ALL_SKY_OPTICS_REFERENCE` | `validation/reference/ecrad/ecckd_32x64_all_sky_tropical_column.nc` |
| `RH_ALL_SKY_OPTICS_PROPERTIES` | `validation/external/ecrad/test/ifs/radiative_properties.nc` |
| `RH_ALL_SKY_OPTICS_OUTPUT_BASENAME` | `ecrad_all_sky_optics_gap_32x64` |
