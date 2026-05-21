# Reduced ecCKD Joint Weight/Block Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_clear_sky_tropical_column |
| Target metric | toa_forcing_max_abs |
| Pressure move count | 4 |
| Starting active move count | 405 |
| Ranked candidate count | 22000 |
| Basis count | 8 |
| Group limit | 8 |
| Block sizes | 2 |
| Probe step | 0.03125 |
| Max table log scale | 0.03125 |
| Max weight logit delta | 0.03125 |
| Residual mode | worst_metric |
| Best ridge lambda | 100 |
| Best direction | fitted |
| Base objective | 8.32190990405 |
| Candidate objective | 9.31858911278 |
| Objective reduction | -0.996679208729 |
| Accepted | false |
| Accepted move count | 0 |
| Raw delta norm | 0.11819885056 |
| Table delta norm | 0.00061305972329 |
| Weight delta norm | 0.0656108816837 |

This diagnostic fits grouped active-entry table amplitudes and shortwave
weight logits in one linearized least-squares system, then accepts the
candidate only after a full nonlinear hard-gate evaluation.
