# Reduced ecCKD Topology Slot Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Radius | 2 |
| Candidates | 33 / 33 |
| Pressure moves | 4 |
| Active moves | 2 |
| Base objective | 8.61303768346 |
| Best objective | 11.1847870777 |
| Best improvement | -2.57174939429 |
| Improved | false |
| Best move | g16 -> g18 |

| Removed | Added | Objective | Improvement | Worst case | Worst metric |
|---:|---:|---:|---:|---|---|
| 1 | 2 | 577.508847895 | -568.895810211 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 1 | 3 | 94.5855704473 | -85.9725327638 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 4 | 2 | 41.7878354069 | -33.1747977235 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 4 | 3 | 138.976819959 | -130.363782276 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 4 | 5 | 136.899693406 | -128.286655722 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 4 | 6 | 55.0849764421 | -46.4719387586 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 9 | 7 | 63.9307432433 | -55.3177055599 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 9 | 8 | 53.715287654 | -45.1022499706 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 9 | 11 | 56.1557724886 | -47.5427348051 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 10 | 8 | 364.216229649 | -355.603191965 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 10 | 11 | 297.476987815 | -288.863950132 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 12 | 11 | 88.0620104679 | -79.4489727844 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 13 | 11 | 60.7462685744 | -52.133230891 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 13 | 15 | 50.154856651 | -41.5418189675 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 14 | 15 | 55.3866719963 | -46.7736343129 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 16 | 15 | 54.5572627242 | -45.9442250408 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 16 | 17 | 13.5944986244 | -4.98146094092 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 16 | 18 | 11.1847870777 | -2.57174939429 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 21 | 19 | 18.5840809668 | -9.97104328332 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 21 | 20 | 12.2518471079 | -3.63880942448 | ecckd_clear_sky_tropical_column | surface_forcing_max_abs |
| 21 | 23 | 237.928541474 | -229.315503791 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 22 | 20 | 29.3468129885 | -20.733775305 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 22 | 23 | 173.43011184 | -164.817074157 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 22 | 24 | 292.235479751 | -283.622442067 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 25 | 23 | 20.1847126818 | -11.5716749983 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 25 | 24 | 19.6005773413 | -10.9875396578 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 25 | 26 | 20.6982471142 | -12.0852094307 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 27 | 26 | 28.0463621009 | -19.4333244174 | ecckd_rcemip_style_column_subset | toa_forcing_max_abs |
| 27 | 29 | 151.968206959 | -143.355169276 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 28 | 26 | 87.0827517617 | -78.4697140783 | ecckd_rcemip_style_column_subset | toa_forcing_max_abs |
| 28 | 29 | 137.237978658 | -128.624940975 | ecckd_rcemip_style_column_subset | surface_forcing_max_abs |
| 30 | 29 | 21.1644372722 | -12.5513995887 | ecckd_clear_sky_tropical_column | heating_rate_rmse |
| 31 | 29 | 29.7922334227 | -21.1791957392 | ecckd_clear_sky_tropical_column | heating_rate_rmse |

This diagnostic keeps optimized local parameter slots fixed while swapping nearby official g-points. It tests whether a neighbor topology can inherit the current slot weights/scales/table moves and improve the nonlinear hard objective before any expensive topology reoptimization.
