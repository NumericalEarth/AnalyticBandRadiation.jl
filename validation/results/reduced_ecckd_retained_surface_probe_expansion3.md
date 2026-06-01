# Reduced ecCKD Retained Surface-Probe Expansion 3

Status: **surface_probe_expansion3_rejected**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | surface_probe3_32_small_tiny |
| Best exact objective | 7.19356726474 |
| Best objective reduction | -4.30801960274e-05 |
| Best TOA forcing | 2.15807017942 W m^-2 |
| Best surface forcing | 2.02800401094 W m^-2 |
| Any accepted | false |

This diagnostic repeats the surface-only linearized residual update after
the first two surface-probe expansions have been promoted, accepting only
strict Pareto-safe full-objective updates.

## Configurations

| Label | Candidates | Probe step | Max log scale | Base objective | Base TOA | Base surface | Best objective | Reduction | TOA | Surface | Pareto safe | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| surface_probe3_32_medium_tiny | 32 | 0.000244140625 | 0.001953125 | 7.19352418454 | 2.15805725536 | 2.02792974588 | 7.19356726838 | -4.30838401169e-05 | 2.15807018051 | 2.02800401744 | false | false |
| surface_probe3_32_small_tiny | 32 | 0.0001220703125 | 0.0009765625 | 7.19352418454 | 2.15805725536 | 2.02792974588 | 7.19356726474 | -4.30801960274e-05 | 2.15807017942 | 2.02800401094 | false | false |
