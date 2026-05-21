# ecRad All-Sky Cloud Sweep

This diagnostic sweeps a small set of all-sky cloud optics knobs and rewrites the candidate NetCDF variables with the best trial at the end. It is diagnostic only; the hard source of truth remains `validation/ecrad_accuracy_gate.jl`.

Best trial: `table_scattering_tripleclouds_alpha_cf1_p2_cloudy_region_lw_scattering_tripleclouds_ifs_aerosol`

| Trial | Worst threshold ratio | TOA forcing | Surface forcing | TOA SW CRE max | Surface SW CRE max | Profile SW CRE max | Heating CRE max |
|---|---:|---:|---:|---:|---:|---:|---:|
| `current_ssa1_g085_075_cf05_scale1` | 466.706 | 140.011874694 | 46.2718498694 | 173.153035643 | 71.3559384487 | 174.704621488 | 25.4254511514 |
| `weaker_cloud_scale075` | 325.961 | 86.6983015798 | 97.7883523008 | 119.839462528 | 42.7577951794 | 120.65128522 | 25.6348055313 |
| `weaker_cloud_scale05` | 523.479 | 132.125333681 | 157.043786356 | 96.9494537342 | 108.600732067 | 159.39220619 | 26.0511078622 |
| `stronger_cloud_scale125` | 587.486 | 176.245789026 | 70.440370381 | 209.386949975 | 97.7598908798 | 211.46537435 | 25.3490601268 |
| `cloud_fraction_linear` | 486.148 | 123.994841703 | 145.844429576 | 156.969695951 | 95.0680919908 | 158.28801945 | 24.8949455105 |
| `cloud_fraction_binary_like` | 811.627 | 243.488147688 | 156.870942902 | 291.414097756 | 209.812363563 | 295.295139755 | 32.7481194409 |
| `less_forward_scattering` | 685.17 | 205.55112653 | 141.828946192 | 241.691307016 | 196.859503313 | 244.530843036 | 25.2192326151 |
| `earlier_ssa099_098` | 614.57 | 184.371149307 | 164.678002115 | 133.550882928 | 191.997522614 | 209.329476249 | 50.5274973078 |
| `table_scattering` | 250.433 | 75.1298410654 | 57.816122558 | 87.0437607612 | 83.3743614503 | 105.385795018 | 17.893024067 |
| `table_scattering_ifs_aerosol_table` | 235.614 | 70.6841889948 | 60.7256831278 | 70.0938431833 | 66.6721096298 | 90.7241234686 | 17.7733851553 |
| `table_scattering_ifs_aerosol_table_fixed_rh08` | 235.183 | 70.5550212898 | 60.605557756 | 72.6938781774 | 69.6021332952 | 93.4028423944 | 17.7541620731 |
| `table_scattering_scale2` | 524.743 | 157.422820542 | 143.913129782 | 187.501807159 | 186.172265096 | 243.641403797 | 11.7035686398 |
| `table_scattering_scale3` | 828.434 | 248.53013789 | 230.038556029 | 278.609124506 | 275.402048286 | 363.253917381 | 17.4457014655 |
| `table_scattering_overlap_maximum` | 1971.8 | 571.116327533 | 591.541359811 | 559.202407837 | 565.983120918 | 587.66481781 | 954.065064322 |
| `table_scattering_overlap_cloudy_region_average` | 1971.8 | 571.116327533 | 591.541359811 | 559.202407837 | 565.983120918 | 587.66481781 | 787.028942946 |
| `table_scattering_overlap_cloudy_region_maximum` | 2226.07 | 571.116327533 | 591.541359811 | 559.202407837 | 565.983120918 | 587.66481781 | 1113.06941727 |
| `table_scattering_adding_cloudy_region` | 1077.76 | 311.580546714 | 323.329331458 | 341.65953333 | 365.588466772 | 408.260110285 | 12.8094104754 |
| `table_scattering_matrix_maximum_cloudy_region` | 626.422 | 187.926519022 | 174.68733371 | 207.945272157 | 227.030446264 | 241.664885627 | 15.2623485933 |
| `table_scattering_matrix_alpha_cloudy_region` | 669.647 | 200.893984465 | 189.216577836 | 220.912737599 | 241.55969039 | 256.094543518 | 15.1391519556 |
| `table_scattering_tripleclouds_alpha_cloudy_region` | 507.395 | 152.218503175 | 151.870374632 | 182.297489792 | 194.129509945 | 211.001605257 | 17.1192948375 |
| `table_scattering_tripleclouds_alpha_cf1_cloudy_region` | 133.987 | 18.9340622719 | 40.1960698476 | 15.9868150555 | 18.1846634478 | 18.1846634478 | 13.9907287521 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region` | 130.718 | 17.8101086443 | 39.2152917143 | 28.9522783952 | 29.7093788302 | 33.0800560182 | 13.9662569992 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_lw_overlap` | 116.662 | 15.6290012807 | 34.998745699 | 28.9522783952 | 29.7093788302 | 33.0800560182 | 5.49908015325 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_lw_tripleclouds` | 107.765 | 18.1324096771 | 32.3293976198 | 28.9522783952 | 29.7093788302 | 33.0800560182 | 0.846012522127 |
| `table_scattering_tripleclouds_alpha_cf1_p2_cloudy_region_lw_tripleclouds` | 107.765 | 14.4002750471 | 32.3293976198 | 15.9868150555 | 18.1846634478 | 18.1846634478 | 0.550087104886 |
| `table_scattering_tripleclouds_alpha_cf1_p2_cloudy_region_ifs_aerosol` | 55.1432 | 12.0937365154 | 16.5429555603 | 0.00852997698621 | 0.00624643419815 | 0.00875834214355 | 13.9309330745 |
| `table_scattering_tripleclouds_alpha_cf1_p2_cloudy_region_lw_tripleclouds_ifs_aerosol` | 15.5126 | 4.65378938075 | 0.345125156235 | 0.00852997698621 | 0.00624643419815 | 0.00875834214355 | 0.506940500565 |
| `table_scattering_tripleclouds_alpha_cf1_p2_cloudy_region_lw_scattering_tripleclouds_ifs_aerosol` | 0.387294 | 0.116188256566 | 0.0128536244544 | 0.00852997698621 | 0.00624643419815 | 0.00875834214355 | 0.112223865109 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_ifs_aerosol` | 52.061 | 14.2484935221 | 15.6183036761 | 12.6458147193 | 11.1152875001 | 17.346398506 | 13.9061445165 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_lw_tripleclouds_ifs_aerosol` | 50.7249 | 15.2174822389 | 11.3223413555 | 12.6458147193 | 11.1152875001 | 17.346398506 | 0.866132777048 |
| `table_scattering_tripleclouds_alpha_cf1_p8_cloudy_region_ifs_aerosol` | 85.5952 | 25.6785497353 | 19.552808346 | 22.8488235824 | 22.1061706377 | 28.9660161762 | 13.8650399449 |
| `table_scattering_tripleclouds_alpha_cf1_sw_p8_lw_p4_cloudy_region_lw_tripleclouds` | 107.765 | 29.8131718876 | 32.3293976198 | 38.5116079999 | 37.5785046689 | 44.8923415127 | 1.01371717378 |
| `table_scattering_tripleclouds_alpha_cf1_sw_p12_lw_p4_cloudy_region_lw_tripleclouds` | 117.646 | 35.293890423 | 32.3293976198 | 43.9923265354 | 42.8509402633 | 49.2530934772 | 1.07084222036 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_lw_tripleclouds_sw_ssa108` | 154.848 | 46.4542617104 | 32.3300751814 | 55.1526978227 | 24.353080952 | 55.5220341215 | 4.75699984172 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_lw_tripleclouds_sw_ssa115` | 188.111 | 56.4331711657 | 32.3305750811 | 65.1316072781 | 23.5033626134 | 65.5392741982 | 6.13576568635 |
| `table_scattering_tripleclouds_alpha_cf1_p8_cloudy_region_lw_tripleclouds` | 107.765 | 29.0797342136 | 32.3293976198 | 38.5116079999 | 37.5785046689 | 44.8923415127 | 1.61195689598 |
| `table_scattering_tripleclouds_alpha_cf1_p12_cloudy_region_lw_tripleclouds` | 114.166 | 34.2497897209 | 32.3293976198 | 43.9923265354 | 42.8509402633 | 49.2530934772 | 1.9071241756 |
| `table_scattering_tripleclouds_alpha_cf1_p4_cloudy_region_lw_scattering` | 135.533 | 19.986349712 | 40.6599552608 | 28.9522783952 | 29.7093788302 | 33.0800560182 | 13.9991256846 |
| `table_scattering_tripleclouds_alpha_cf1_p8_cloudy_region` | 125.524 | 28.838737824 | 37.6571352423 | 38.5116079999 | 37.5785046689 | 44.8923415127 | 13.9256471532 |
| `table_scattering_tripleclouds_alpha_p8_cloudy_region` | 551.212 | 165.359591803 | 165.363498177 | 195.438578419 | 207.62263349 | 227.499847411 | 16.8624489972 |
| `table_scattering_tripleclouds_alpha_p12_cloudy_region` | 561.936 | 168.434339164 | 168.580681509 | 198.51332578 | 210.839816822 | 231.186547993 | 16.7405308205 |

## Best Configuration

| Environment variable | Value |
|---|---|
| `RH_AEROSOL_DELTA_EDDINGTON_SCALE` | `true` |
| `RH_AEROSOL_OPTICS` | `true` |
| `RH_CANDIDATE_GAS_OPTICS` | `official_ecckd` |
| `RH_CLOUD_EFFECTIVE_RADIUS_OPTICS` | `true` |
| `RH_CLOUD_FRACTION_EXPONENT` | `1.0` |
| `RH_CLOUD_INHOM_OVERLAP_EXPONENT` | `2.0` |
| `RH_CLOUD_LW_TABLE_EXTINCTION_AS_ABSORPTION` | `false` |
| `RH_CLOUD_OVERLAP_LONGWAVE` | `true` |
| `RH_CLOUD_OVERLAP_LONGWAVE_RULE` | `tripleclouds_alpha` |
| `RH_CLOUD_OVERLAP_RULE` | `tripleclouds_alpha` |
| `RH_CLOUD_OVERLAP_SHORTWAVE` | `true` |
| `RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS` | `true` |
| `RH_CLOUD_SCATTERING_DELTA_EDDINGTON_AVERAGE` | `true` |
| `RH_CLOUD_SCATTERING_DELTA_EDDINGTON_SCALE` | `true` |
| `RH_CLOUD_SCATTERING_LAYER_EFFECTIVE_RADIUS` | `true` |
| `RH_CLOUD_SCATTERING_MAPPING_METHOD` | `ecrad` |
| `RH_CLOUD_SCATTERING_SW_SSA_SCALE` | `1.0` |
| `RH_CLOUD_SCATTERING_TABLE_OPTICS` | `true` |
| `RH_CLOUD_SCATTERING_THICK_AVERAGING` | `true` |
| `RH_ICE_CLOUD_EFFECTIVE_RADIUS_SW_SCALE` | `1.0` |
| `RH_ICE_CLOUD_LW_MASS_ABSORPTION` | `50` |
| `RH_ICE_CLOUD_SCATTERING_SW_SCALE` | `1.0` |
| `RH_ICE_CLOUD_SW_SCATTERING_ASYMMETRY` | `0.75` |
| `RH_ICE_CLOUD_SW_SINGLE_SCATTERING_ALBEDO` | `1.0` |
| `RH_IFS_AEROSOL_LAYER_RELATIVE_HUMIDITY` | `true` |
| `RH_IFS_AEROSOL_RELATIVE_HUMIDITY` | `0.8` |
| `RH_IFS_AEROSOL_TABLE_OPTICS` | `true` |
| `RH_LIQUID_CLOUD_EFFECTIVE_RADIUS_SW_SCALE` | `1.0` |
| `RH_LIQUID_CLOUD_LW_MASS_ABSORPTION` | `100` |
| `RH_LIQUID_CLOUD_SCATTERING_SW_SCALE` | `1.0` |
| `RH_LIQUID_CLOUD_SW_SCATTERING_ASYMMETRY` | `0.85` |
| `RH_LIQUID_CLOUD_SW_SINGLE_SCATTERING_ALBEDO` | `1.0` |
| `RH_LW_CLOUD_SCATTERING` | `true` |
