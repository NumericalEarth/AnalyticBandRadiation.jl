# Official/Reduced ecCKD Gas-Optics Training

Status: **passed**

| Field | Value |
|---|---:|
| Trainable shortwave g-points | 16 |
| Parameter count | 48 |
| Iterations | 2 |
| Initial objective | 214.264529894 |
| Final objective | 213.306642558 |
| Objective reduction | 0.957887335617 |
| Final objective / initial objective | 0.995529417136 |
| Final objective / hard target | 213.306642558 |
| Hard accuracy target met | false |
| Reactant check | passed |
| Enzyme check | passed |

This is the official/reduced ecCKD training-path artifact: it demonstrates objective construction, trainable parameters, Reactant/Enzyme checks, and deterministic objective reduction on official ecCKD references. It does not close the reduced-accuracy gate until final_objective_target_ratio <= 1.

Next required work: Replace the current per-g-point scaling search with a richer reduced k-distribution coefficient/topology optimization; the current 48-parameter search remains far above the hard-gate target.
