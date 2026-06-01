# Reduced ecCKD Exact Weight Refit

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Iterations requested | 8 |
| Iterations completed | 1 |
| Initial objective | 8.32190990405 |
| Final objective | 8.32190990405 |
| Objective reduction | 0 |
| Accepted | false |
| Initial worst TOA forcing | 2.49657297122 W m^-2 |
| Initial worst surface forcing | 2.38029831045 W m^-2 |
| Final worst TOA forcing | 2.49657297122 W m^-2 |
| Final worst surface forcing | 2.38029831045 W m^-2 |

This diagnostic refits only the 16 shortwave weights by coordinate
searching softmax logits against the exact clean ecCKD hard objective.
It is accepted only if the exact nonlinear model improves.
