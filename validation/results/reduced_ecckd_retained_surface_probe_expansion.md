# Reduced ecCKD Retained Surface-Probe Expansion

Status: **surface_probe_expansion_rejected**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | surface_probe_32_small_tiny |
| Best exact objective | 7.18884381786 |
| Best objective reduction | -5.74426125395e-05 |
| Best TOA forcing | 2.15665314536 W m^-2 |
| Best surface forcing | 2.01228897428 W m^-2 |
| Any accepted | false |

This diagnostic uses the current retained base with a surface-only
linearized residual and accepts only strict Pareto-safe full-objective
updates that do not regress either TOA or surface forcing.

## Configurations

| Label | Candidates | Probe step | Max log scale | Base objective | Base TOA | Base surface | Best objective | Reduction | TOA | Surface | Pareto safe | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| surface_probe_32_medium_tiny | 32 | 0.000244140625 | 0.001953125 | 7.18878637525 | 2.15663591257 | 2.01218723328 | 7.18884382221 | -5.7446958408e-05 | 2.15665314666 | 2.01228898204 | false | false |
| surface_probe_32_small_tiny | 32 | 0.0001220703125 | 0.0009765625 | 7.18878637525 | 2.15663591257 | 2.01218723328 | 7.18884381786 | -5.74426125395e-05 | 2.15665314536 | 2.01228897428 | false | false |
