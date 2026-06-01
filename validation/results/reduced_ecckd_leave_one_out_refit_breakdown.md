# Reduced ecCKD Leave-One-Out Refit Breakdown

Status: **failed_threshold**

- Omitted SW g-point: 25
- Saved refit objective: 12.332454488
- Recomputed objective: 12.332454488
- Worst case: `ecckd_clear_sky_tropical_column`
- Limiting metric: `heating_rate_rmse` = 0.616622724398 / 0.05 = 12.332454488×

## Compared Leave-One-Out Rows

| Omitted SW g-point | Objective | Worst case | Limiting metric | Value | Threshold | Ratio |
|---:|---:|---|---|---:|---:|---:|
| 25 | 12.332454488 | `ecckd_clear_sky_tropical_column` | `heating_rate_rmse` | 0.616622724398 | 0.05 | 12.332454488 |
| 23 | 1.62442314379 | `ecckd_clear_sky_tropical_column` | `heating_rate_rmse` | 0.0812211571893 | 0.05 | 1.62442314379 |

## Metric Ratios

| Case | Rank | Metric | Value | Threshold | Ratio |
|---|---:|---|---:|---:|---:|
| `ecckd_clear_sky_tropical_column` | 1 | `heating_rate_rmse` | 0.616622724398 | 0.05 | 12.332454488 |
| `ecckd_clear_sky_tropical_column` | 2 | `heating_rate_max_abs` | 6.10463644722 | 0.5 | 12.2092728944 |
| `ecckd_clear_sky_tropical_column` | 3 | `surface_forcing` | 0.0209399713759 | 0.3 | 0.0697999045864 |
| `ecckd_clear_sky_tropical_column` | 4 | `toa_forcing` | 0.0138484135641 | 0.3 | 0.0461613785471 |
| `ecckd_clear_sky_tropical_column` | 5 | `sw_down_rmse` | 0.0331309671129 | 1 | 0.0331309671129 |
| `ecckd_clear_sky_tropical_column` | 6 | `sw_down_max_abs` | 0.0455293915809 | 5 | 0.00910587831618 |
| `ecckd_clear_sky_tropical_column` | 7 | `sw_up_rmse` | 0.0087839735313 | 1 | 0.0087839735313 |
| `ecckd_clear_sky_tropical_column` | 8 | `sw_up_max_abs` | 0.0137202363576 | 5 | 0.00274404727152 |
| `ecckd_rcemip_style_column_subset` | 1 | `heating_rate_max_abs` | 6.10463644722 | 0.5 | 12.2092728944 |
| `ecckd_rcemip_style_column_subset` | 2 | `heating_rate_rmse` | 0.500562631018 | 0.05 | 10.0112526204 |
| `ecckd_rcemip_style_column_subset` | 3 | `surface_forcing` | 0.0242811668201 | 0.3 | 0.0809372227335 |
| `ecckd_rcemip_style_column_subset` | 4 | `toa_forcing` | 0.0207035096073 | 0.3 | 0.069011698691 |
| `ecckd_rcemip_style_column_subset` | 5 | `sw_down_rmse` | 0.0268112083695 | 1 | 0.0268112083695 |
| `ecckd_rcemip_style_column_subset` | 6 | `sw_up_rmse` | 0.00930730931708 | 1 | 0.00930730931708 |
| `ecckd_rcemip_style_column_subset` | 7 | `sw_down_max_abs` | 0.0455293915809 | 5 | 0.00910587831618 |
| `ecckd_rcemip_style_column_subset` | 8 | `sw_up_max_abs` | 0.0214792939209 | 5 | 0.00429585878419 |
