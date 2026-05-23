# ecCKD Published Recovery Target

Status: **published_recovery_target_ready**

Optimizer-only delta rule: A recovery run may vary optimizer settings, schedules, and initialization seeds, but must keep the published CKDMIP/ecCKD inputs, objective terms, trainable arrays, and evaluation cases fixed.

## Acceptance Metrics

| Metric | Threshold |
|---|---:|
| `final_objective_target_ratio_max` | 1.05 |
| `forcing_error_regression_margin_w_m2` | 0.03 |
| `heating_rmse_regression_margin_k_day` | 0.005 |
| `optical_depth_log_rmse_max` | 0.02 |
| `weight_l1_relative_error_max` | 0.02 |

## Published Targets

| File | Kind | Bands | G-points | Coefficient arrays | Coefficient parameters | Support arrays | Support parameters |
|---|---|---:|---:|---:|---:|---:|---:|
| `ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc` | longwave | 1 | 32 | 8 | 193344 | 3 | 17856 |
| `ecckd-1.0_sw_climate_rgb-32b_ckd-definition.nc` | shortwave | 5 | 32 | 6 | 172992 | 4 | 31936 |
| `ecckd-1.2_lw_climate_narrow-64b_ckd-definition.nc` | longwave | 13 | 64 | 8 | 386688 | 3 | 35712 |
| `ecckd-1.2_sw_climate_window-64b_ckd-definition.nc` | shortwave | 19 | 64 | 6 | 345984 | 4 | 63872 |
| `ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc` | shortwave | 5 | 32 | 6 | 172992 | 4 | 31936 |
| `ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc` | shortwave | 44 | 96 | 6 | 518976 | 4 | 95808 |

## Primary 32-G Targets

- `ecckd-1.0_lw_climate_fsck-32b_ckd-definition.nc`: longwave, 1 bands, 32 g-points, 193344 coefficient parameters.
- `ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc`: shortwave, 5 bands, 32 g-points, 172992 coefficient parameters.

## Next Required Work

- Choose one primary target, preferably the 32-g shortwave RGB model, and initialize the trainable arrays from the same parameterization used by ecCKD.
- Run the original-objective optimizer with only optimizer settings varied.
- Write the recovered CKD-definition and evaluate it with ecckd_recovery_metrics plus the original-objective flux/heating criteria.
