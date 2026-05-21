# ecRad Accuracy Diagnostics

Gate status: **passed**

This report ranks current `radiative_heating_*` errors by threshold exceedance. It is diagnostic only; the hard pass/fail source of truth remains `validation/ecrad_accuracy_gate.jl`.

Failed metric count: 0

## Case Summary

| Case | Failed metrics | Worst threshold ratio |
|---|---:|---:|
| ecckd_clear_sky_tropical_column | 0 | 0.182542 |
| ecckd_all_sky_tropical_column | 0 | 0.387294 |
| ecckd_rcemip_style_column_subset | 0 | 0.182542 |

## Worst Metrics

| Case | Metric | Value | Threshold | Ratio |
|---|---|---:|---:|---:|
| ecckd_all_sky_tropical_column | `toa_forcing_abs_error` | 0.116188256566 W m^-2 | 0.3 W m^-2 | 0.387294 |
| ecckd_clear_sky_tropical_column | `heating_rate_max_abs` | 0.0912709758456 K day^-1 | 0.5 K day^-1 | 0.182542 |
| ecckd_rcemip_style_column_subset | `heating_rate_max_abs` | 0.0912709758456 K day^-1 | 0.5 K day^-1 | 0.182542 |
| ecckd_all_sky_tropical_column | `heating_rate_max_abs` | 0.0768932791495 K day^-1 | 0.5 K day^-1 | 0.153787 |
| ecckd_all_sky_tropical_column | `heating_rate_rmse` | 0.00575848422016 K day^-1 | 0.05 K day^-1 | 0.11517 |
| ecckd_clear_sky_tropical_column | `heating_rate_rmse` | 0.00575159066875 K day^-1 | 0.05 K day^-1 | 0.115032 |
| ecckd_rcemip_style_column_subset | `heating_rate_rmse` | 0.0045715395775 K day^-1 | 0.05 K day^-1 | 0.0914308 |
| ecckd_clear_sky_tropical_column | `surface_forcing_abs_error` | 0.0140396962501 W m^-2 | 0.3 W m^-2 | 0.046799 |
| ecckd_rcemip_style_column_subset | `surface_forcing_abs_error` | 0.0140396962501 W m^-2 | 0.3 W m^-2 | 0.046799 |
| ecckd_all_sky_tropical_column | `surface_forcing_abs_error` | 0.0128536244544 W m^-2 | 0.3 W m^-2 | 0.0428454 |
| ecckd_all_sky_tropical_column | `lw_up_rmse` | 0.0375781506742 W m^-2 | 1 W m^-2 | 0.0375782 |
| ecckd_all_sky_tropical_column | `lw_up_max_abs` | 0.147512483745 W m^-2 | 5 W m^-2 | 0.0295025 |
