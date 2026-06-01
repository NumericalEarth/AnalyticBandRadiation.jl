# Reduced ecCKD Global Block Linearized Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_clear_sky_tropical_column |
| Target metric | toa_forcing_max_abs |
| Pressure move count | 4 |
| Starting active move count | 405 |
| Ranked candidate count | 22000 |
| Basis count | 48 |
| Group limit | 24 |
| Block sizes | 2, 4 |
| Probe step | 0.0625 |
| Max log scale | 0.03125 |
| Best ridge lambda | 0.01 |
| Base objective | 8.39187324383 |
| Candidate objective | 8.39467524724 |
| Objective reduction | -0.00280200340474 |
| Accepted | false |
| Accepted move count | 0 |
| Raw delta norm | 166.72317203 |
| Clamped delta norm | 0.154716087081 |

This diagnostic fits a joint linearized update over grouped active-entry
table bases. It differs from scalar block refinement by solving for all
basis amplitudes at once, then accepting the candidate only after a full
nonlinear hard-gate evaluation.
