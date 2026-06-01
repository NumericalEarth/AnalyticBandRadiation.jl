# Reduced ecCKD Retained Structural Pareto Probe

Status: **pareto_probe_rejected**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | current_residual_tiny_no_rayleigh |
| Best exact objective | 7.19352418454 |
| Best objective reduction | 0 |
| Best TOA forcing | 2.15805725536 W m^-2 |
| Best surface forcing | 2.02838387132 W m^-2 |
| Any accepted | false |
| Any Pareto-safe accepted | false |

This diagnostic probes the current retained structural base after the
four retained structural continuations. It records whether another tiny
constrained table step is available without repeating the retained
one-off continuation chain.

## Configurations

| Label | Scope | Rayleigh | Base objective | Best objective | Reduction | TOA | Surface | TOA regressed | Surface regressed | Accepted |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| current_exact_objective_tiny_rayleigh | all_global_objective_probe | true | 7.19352418454 | 7.19485514125 | -0.00133095670468 | 2.15845654237 | 2.03049539621 | true | true | false |
| current_residual_tiny_no_rayleigh | all_global_residual_probe | false | 7.19352418454 | 7.19352418454 | 0 | 2.15805725536 | 2.02838387132 | false | true | false |
