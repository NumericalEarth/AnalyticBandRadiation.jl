# Reduced ecCKD Gas Pressure-Band Linearized Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_clear_sky_tropical_column |
| Target metric | toa_forcing_max_abs |
| Target boundary | toa |
| Target column | 7 |
| Target residual | 2.49621126044 W m^-2 |
| Fitted candidates | 24 |
| Probe step | 0.00048828125 |
| Max log scale | 0.00048828125 |
| Best ridge lambda | 1 |
| Best direction | fitted |
| Base objective | 8.32070420146 |
| Candidate objective | 8.32072913222 |
| Objective reduction | -2.49307637059e-05 |
| Accepted | false |
| Accepted move count | 0 |
| Raw delta norm | 6.9992888185 |
| Clamped delta norm | 0.00202644778949 |

This diagnostic solves a ridge-regularized linearized update over
gas-specific pressure-band scales for the current worst boundary
residual, then accepts only after full nonlinear hard-gate evaluation.
