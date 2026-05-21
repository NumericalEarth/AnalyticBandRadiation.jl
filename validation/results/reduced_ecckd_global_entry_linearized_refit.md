# Reduced ecCKD Global Entry Linearized Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Target case | ecckd_rcemip_style_column_subset |
| Target metric | toa_forcing_max_abs |
| Pressure move count | 4 |
| Starting active move count | 8 |
| Ranked candidate count | 34048 |
| Fitted candidate count | 24 |
| Probe step | 0.0625 |
| Max log scale | 0.25 |
| Best ridge lambda | 1e-06 |
| Base objective | 8.57260015149 |
| Candidate objective | 8.57260015149 |
| Objective reduction | 0 |
| Accepted | false |
| Accepted move count | 0 |
| Raw delta norm | 7158.12256329 |
| Clamped delta norm | 1.14564392374 |

This diagnostic starts from the best available pressure/table-active
state and fits a joint linearized update over the highest-priority
all-slot active table entries. The update is accepted only if the full
hard-gate objective improves after evaluating the nonlinear model.
