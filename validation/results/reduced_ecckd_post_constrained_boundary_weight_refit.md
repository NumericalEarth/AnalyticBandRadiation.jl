# Reduced ecCKD Post-Constrained Boundary Weight Refit

Status: **boundary_weight_refit_improved**

| Field | Value |
|---|---:|
| Iterations requested | 20 |
| Iterations completed | 2 |
| Initial boundary objective | 7.2974317537 |
| Final boundary objective | 7.29742558538 |
| Boundary objective reduction | 6.16832143852e-06 |
| Initial full objective | 7.2974317537 |
| Final full objective | 7.29742558538 |
| Accepted | true |
| Initial worst TOA forcing | 2.18922952611 W m^-2 |
| Initial worst surface forcing | 2.18919927634 W m^-2 |
| Final worst TOA forcing | 2.18922717474 W m^-2 |
| Final worst surface forcing | 2.18922767561 W m^-2 |

This diagnostic refits only shortwave quadrature weights against the
boundary-forcing max objective after cumulative constrained table moves.
It is diagnostic and is not the canonical reduced model unless all hard
thresholds pass.
