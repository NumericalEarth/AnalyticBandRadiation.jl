# Reduced ecCKD Targeted Entry Linearized Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Pressure move count | 4 |
| Ranked candidate count | 1155 |
| Fitted candidate count | 64 |
| Probe step | 0.0625 |
| Max log scale | 0.25 |
| Best ridge lambda | 1e-06 |
| Base objective | 8.61316394314 |
| Candidate objective | 8.61303778873 |
| Objective reduction | 0.00012615441913 |
| Accepted | true |
| Accepted move count | 64 |
| Raw delta norm | 7676.97230482 |
| Clamped delta norm | 2 |

This diagnostic fits a joint linearized update over the highest-priority
active table entries. The update is accepted only if the full hard-gate
objective improves after evaluating the nonlinear model.
