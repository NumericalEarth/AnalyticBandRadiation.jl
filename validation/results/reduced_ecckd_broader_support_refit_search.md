# Reduced ecCKD Broader Support-Plus-Refit Search

Status: **support_refit_rejected**

| Field | Value |
|---|---:|
| Objective target | 1 |
| Current objective | 7.00640773953 |
| Best objective | 45.0350379919 |
| Best objective reduction | -38.0286302524 |
| Best passed hard objective | false |
| Radius | 2 |
| Prefilter candidates | 8 / 33 |
| Candidates evaluated | 2 / 33 |
| Pressure bands | 4 |
| Pressure partition | log_pressure |
| Include Rayleigh | false |
| Best support move | g4 -> g2 |

This artifact is intentionally stricter than prior support-only scans:
each alternate support is evaluated after replaying the current composed
retained table/component/gas-pressure chain and then running one bounded
static-gas/H2O pressure-component refit against the full reduced hard
objective.

## Candidate Rows

| Move | Base objective | Refit best objective | Selected objective | TOA | Surface | Heating RMSE | Accepted moves |
|---|---:|---:|---:|---:|---:|---:|---:|
| g4 -> g2 | 45.0350379919 | 44.7942188596 | 45.0350379919 | 5.41062119486 | 13.5105113976 | 0.354548437281 | 0 |
| g4 -> g6 | 50.7860997415 | 50.4799095623 | 50.7860997415 | 3.90538360112 | 15.2358299224 | 0.399675465158 | 0 |
