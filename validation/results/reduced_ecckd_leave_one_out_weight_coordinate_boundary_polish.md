# Reduced ecCKD Leave-One-Out Weight Coordinate Boundary Polish

Status: **weight_coordinate_boundary_polish_passed**

| Field | Value |
|---|---:|
| Omitted SW g-point | 23 |
| Model | 32x31 |
| Coordinate count | 31 |
| Max iterations | 20 |
| Completed iterations | 17 |
| Accepted moves | 17 |
| Initial objective | 1.08338530362 |
| Final objective | 0.999422334049 |
| Objective reduction | 0.0839629695705 |
| Final TOA forcing | 0.299826700215 W m^-2 |
| Final surface forcing | 0.278191612233 W m^-2 |
| Final boundary forcing | 0.299826700215 W m^-2 |
| Final heating RMSE | 0.0498769910608 K day^-1 |

This diagnostic starts from the boundary-tight continuation row and
uses very small normalized quadrature-weight moves. It stops once
the exact hard objective reaches the pass threshold or no accepted
boundary-safe descent move remains.

## Iterations

| Iteration | Accepted | G-point | Log scale | Objective | Reduction | Boundary forcing | Heating RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | true | 25 | 0.001 | 1.07738459671 | 0.00600070690824 | 0.299456385968 | 0.0538692298356 |
| 2 | true | 25 | 0.001 | 1.07149146795 | 0.00589312876457 | 0.299479357441 | 0.0535745733973 |
| 3 | true | 25 | 0.001 | 1.06570814018 | 0.00578332777107 | 0.299502351896 | 0.0532854070088 |
| 4 | true | 25 | 0.001 | 1.06003685235 | 0.00567128782845 | 0.299525369356 | 0.0530018426174 |
| 5 | true | 25 | 0.001 | 1.0544798573 | 0.00555699504668 | 0.299548409842 | 0.052723992865 |
| 6 | true | 25 | 0.001 | 1.0490394198 | 0.0054404375015 | 0.299571473379 | 0.0524519709899 |
| 7 | true | 25 | 0.001 | 1.04371781433 | 0.00532160546973 | 0.299594559989 | 0.0521858907165 |
| 8 | true | 25 | 0.001 | 1.03851732276 | 0.00520049156755 | 0.299617669696 | 0.0519258661381 |
| 9 | true | 25 | 0.001 | 1.03344023207 | 0.00507709069203 | 0.299640802521 | 0.0516720116035 |
| 10 | true | 25 | 0.001 | 1.02848883138 | 0.004951400689 | 0.29966395849 | 0.051424441569 |
| 11 | true | 25 | 0.001 | 1.02366540963 | 0.00482342175112 | 0.299687137625 | 0.0511832704815 |
| 12 | true | 25 | 0.001 | 1.01897225237 | 0.00469315725465 | 0.299710339949 | 0.0509486126187 |
| 13 | true | 25 | 0.001 | 1.01441163915 | 0.00456061322101 | 0.299733565484 | 0.0507205819577 |
| 14 | true | 25 | 0.001 | 1.00998584007 | 0.00442579908466 | 0.299756814256 | 0.0504992920035 |
| 15 | true | 25 | 0.001 | 1.00569711281 | 0.00428872726181 | 0.299780086286 | 0.0502848556404 |
| 16 | true | 25 | 0.001 | 1.00154769907 | 0.00414941373928 | 0.299803381598 | 0.0500773849534 |
| 17 | true | 25 | 0.001 | 0.999422334049 | 0.00212536501913 | 0.299826700215 | 0.0498769910608 |
