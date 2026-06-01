# Reduced ecCKD Nonlocal Support-Plus-Refit Search

Status: **nonlocal_support_refit_rejected**

| Field | Value |
|---|---:|
| Objective target | 1 |
| Current objective | 7.00640773953 |
| Best objective | 45.0350379919 |
| Best objective reduction | -38.0286302524 |
| Best passed hard objective | false |
| Candidate sources | 4 |
| Prefilter candidates | 4 / 4 |
| Candidates evaluated | 4 / 4 |
| Pressure bands | 4 |
| Pressure partition | log_pressure |
| Include Rayleigh | false |
| Best label | hardgate_subset_best |

This artifact tests nonlocal 16-g supports already identified by the
random, swap, continuation, and hardgate subset searches. Each support
is evaluated after replaying the current composed table/component/
gas-pressure chain and then running one bounded static-gas/H2O
pressure-component refit against the full reduced hard objective.

## Candidate Rows

| Label | Base objective | Refit best objective | Selected objective | TOA | Surface | Heating RMSE | Accepted moves |
|---|---:|---:|---:|---:|---:|---:|---:|
| hardgate_subset_best | 45.0350379919 | 44.7942188596 | 45.0350379919 | 5.41062119486 | 13.5105113976 | 0.354548437281 | 0 |
| support_swap_best | 1375.13330949 | 1374.38822497 | 1375.13330949 | 183.280058873 | 412.539992846 | 41.5707631954 | 0 |
| support_swap_continuation_best | 70366.6855162 | 70279.8096622 | 70366.6855162 | 183.280143057 | 412.613762653 | 3518.33427581 | 0 |
| random_support_best | 79316.6222921 | 79219.020283 | 79316.6222921 | 314.434117001 | 422.004875708 | 3965.83111461 | 0 |
