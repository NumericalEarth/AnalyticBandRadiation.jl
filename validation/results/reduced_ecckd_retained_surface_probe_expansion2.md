# Reduced ecCKD Retained Surface-Probe Expansion 2

Status: **surface_probe_expansion2_rejected**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | surface_probe2_32_small_tiny |
| Best exact objective | 7.19271692514 |
| Best objective reduction | -4.5715205867e-05 |
| Best TOA forcing | 2.15781507754 W m^-2 |
| Best surface forcing | 2.02376849152 W m^-2 |
| Any accepted | false |

This diagnostic repeats the surface-only linearized residual update after
the first surface-probe expansion has been promoted, accepting only
strict Pareto-safe full-objective updates.

## Configurations

| Label | Candidates | Probe step | Max log scale | Base objective | Base TOA | Base surface | Best objective | Reduction | TOA | Surface | Pareto safe | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| surface_probe2_32_medium_tiny | 32 | 0.000244140625 | 0.001953125 | 7.19267120994 | 2.15780136298 | 2.02368917152 | 7.19271692892 | -4.5718978896e-05 | 2.15781507867 | 2.02376849816 | false | false |
| surface_probe2_32_small_tiny | 32 | 0.0001220703125 | 0.0009765625 | 7.19267120994 | 2.15780136298 | 2.02368917152 | 7.19271692514 | -4.5715205867e-05 | 2.15781507754 | 2.02376849152 | false | false |
