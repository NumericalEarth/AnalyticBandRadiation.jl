# Official/Reduced ecCKD Gas-Optics Training

Status: **partial**

| Field | Value |
|---|---:|
| Trainable shortwave g-points | 16 |
| Parameter count | 48 |
| Iterations | 27 |
| Initial objective | 214.264529894 |
| Final objective | 8.60500399071 |
| Objective reduction | 205.659525903 |
| Final objective / initial objective | 0.0401606555923 |
| Final objective / hard target | 8.60500399071 |
| Hard accuracy target met | false |
| Recovery status | optimizer_improved_but_target_not_met |
| Reactant check | passed |
| Enzyme check | passed |

This is the official/reduced ecCKD training-path artifact: it demonstrates objective construction, trainable parameters, Reactant/Enzyme checks, and deterministic objective reduction on official ecCKD references. Status is partial until final_objective_target_ratio <= 1 and a published model is recovered quantitatively.

Next required work: Move beyond bounded pressure-band table scales: run a stronger joint coefficient/table optimizer against flux and heating residuals or jointly optimize the reduced quadrature definition; the current table-refined 48-parameter path remains far above the hard-gate target.
