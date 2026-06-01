# Reduced ecCKD Retained Topology Constrained Optimizer

Status: **constrained_table_optimizer_rejected**

| Field | Value |
|---|---:|
| Base mode | retained_topology |
| Candidate scope | all_global_residual_probe |
| Residual mode | boundary |
| SW indices | 1, 4, 9, 10, 12, 13, 14, 16, 19, 22, 26, 27, 28, 30, 31, 32 |
| Candidate count | 16 |
| Probe step | 0.00390625 |
| Max log scale | 0.03125 |
| Base objective | 8.80904002845 |
| Best exact objective | 9.49940088339 |
| Best TOA forcing error | 2.27359978558 W m^-2 |
| Best surface forcing error | 2.84982026502 W m^-2 |
| Pareto safe | false |
| Accepted | false |
| Proposed moves | 16 |

This diagnostic checks whether the best rejected retained-topology pair can
be rescued by a small local constrained table-entry optimization before
considering broader topology or coefficient-table parameterizations.
