# Reduced ecCKD Retained Objective-Probe Expansion

Status: **objective_probe_expansion_rejected**

| Field | Value |
|---|---:|
| Configurations | 2 |
| Best label | objective_probe_24_medium_tiny |
| Best exact objective | 7.14795534014 |
| Best objective reduction | 0.00535018525099 |
| Best TOA forcing | 2.14438660204 W m^-2 |
| Best surface forcing | 2.05760969562 W m^-2 |
| Any accepted | false |

This diagnostic expands the retained constrained-table objective probe
beyond the 12-candidate tiny Pareto probe. It uses exact hard-objective
single-coordinate probes to select 24 table/Rayleigh candidates, then
solves the bounded linearized residual system and exact-evaluates the
candidate with strict Pareto-safe full-objective acceptance.

## Configurations

| Label | Candidates | Probe step | Max log scale | Base objective | Base TOA | Base surface | Best objective | Reduction | TOA | Surface | Pareto safe | Accepted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| objective_probe_24_medium_tiny | 24 | 0.000244140625 | 0.001953125 | 7.15330552539 | 2.14599165762 | 2.03415646401 | 7.14795534014 | 0.00535018525099 | 2.14438660204 | 2.05760969562 | false | false |
| objective_probe_24_small_tiny | 24 | 0.0001220703125 | 0.0009765625 | 7.15330552539 | 2.14599165762 | 2.03415646401 | 7.1506284651 | 0.00267706028571 | 2.14518853953 | 2.04586875075 | false | false |
