# Reduced ecCKD Leave-One-Out Weight Coordinate Descent

Status: **weight_coordinate_descent_improved**

| Field | Value |
|---|---:|
| Omitted SW g-point | 23 |
| Model | 32x31 |
| Max iterations | 6 |
| Completed iterations | 6 |
| Accepted moves | 6 |
| Initial objective | 1.62442314379 |
| Final objective | 1.26634555306 |
| Objective reduction | 0.358077590727 |
| Final TOA forcing | 0.299874474241 W m^-2 |
| Final surface forcing | 0.278151667767 W m^-2 |
| Final boundary forcing | 0.299874474241 W m^-2 |
| Final heating RMSE | 0.063317277653 K day^-1 |

This diagnostic greedily repeats the exact one-coordinate normalized
weight move that improves the hard objective, keeps boundary forcing
under the 0.3 W m^-2 cap, and does not regress heating RMSE.

## Iterations

| Iteration | Accepted | G-point | Log scale | Objective | Reduction | Boundary forcing | Heating RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | true | 28 | -0.01 | 1.48556457006 | 0.138858573726 | 0.254747858736 | 0.074278228503 |
| 2 | true | 11 | -0.0316227766017 | 1.38732357452 | 0.0982409955446 | 0.218255120345 | 0.0693661787258 |
| 3 | true | 8 | -0.0316227766017 | 1.33252528987 | 0.0547982846466 | 0.271310219162 | 0.0666262644935 |
| 4 | true | 29 | -0.01 | 1.28978350614 | 0.0427417837294 | 0.295505033736 | 0.064489175307 |
| 5 | true | 1 | -0.01 | 1.27259497824 | 0.0171885279016 | 0.297702472506 | 0.0636297489119 |
| 6 | true | 11 | -0.00316227766017 | 1.26634555306 | 0.00624942517928 | 0.299874474241 | 0.063317277653 |
