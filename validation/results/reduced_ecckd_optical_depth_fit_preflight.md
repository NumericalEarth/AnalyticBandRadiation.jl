# Reduced ecCKD Optical-Depth Fit Preflight

Status: **optical_depth_refit_target_ready**

| Field | Value |
|---|---:|
| Shortwave g-points | 16 |
| Optical-depth samples | 92064 |
| Baseline optical-depth RMSE | 0.366488363887 |
| Fitted optical-depth RMSE | 0.197241557253 |
| Component-fitted optical-depth RMSE | 0.165192904749 |
| Relative RMSE reduction | 0.461807 |
| Component relative RMSE reduction | 0.549255 |
| Baseline relative RMSE | 1.54615186044 |
| Fitted relative RMSE | 0.832128467786 |
| Component fitted relative RMSE | 0.69692067246 |
| Scale range | 0.0143597 to 100 |
| Component scale range | 0 to 100 |
| Baseline flux objective | 214.264529894 |
| Scaled flux objective | 1062.45269062 |
| Scaled flux objective improved | false |
| Component-scaled flux objective | 1182.45978193 |
| Component-scaled flux objective improved | false |
| Coefficient-table fit parameters per g | 6361 |
| Physical projected table optical-depth RMSE | 3.18762096555e-17 |
| Physical projected table flux objective | 549.944116676 |
| Coefficient-table raw LS optical-depth RMSE | 0.0035778590814 |
| Coefficient-table clipped-model optical-depth RMSE | 2.72779995305 |
| Coefficient-table fit clipped parameters | 13972 |
| Coefficient-table fit flux objective | 1395.50127989 |
| Coefficient-table fit flux objective improved | false |

This artifact is an optical-depth training preflight, not a flux acceptance result. It compares the current weighted-greedy 16-g subset to a weighted 32-to-16 cumulative-k projection of the official shortwave ecCKD optical depths.

The fitted per-g optical-depth scales are diagnostic targets only. The coefficient-table fit is the first table-level refit path: it solves for shortwave absorption, dynamic H2O, and Rayleigh table entries against the same optical-depth target, records the raw unconstrained fit separately from the clipped physical table, and then evaluates the resulting fluxes.

Next required work: The first real coefficient-table least-squares refit is wired, but unconstrained coefficients require heavy nonnegative clipping and these optical-depth targets are still not enough for the flux hard gate; next use a constrained table optimizer against flux/heating residuals or optimize the reduced quadrature definition jointly with the table entries.
