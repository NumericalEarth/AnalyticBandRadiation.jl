# Reduced ecCKD Optimization Preflight

Status: **preflight_ready**

| Field | Value |
|---|---:|
| Trainable shortwave g-points | 16 |
| Parameter count | 48 |
| Initial normalized objective | 214.264529894 |
| Directional derivative | -47.9336647895 |
| Best trial objective | 213.785388667 |
| Best trial step | 0.01 |
| Optimization smoke final objective | 213.785388667 |
| Optimization smoke reduction | 0.47914122634 |
| Directional search final objective | 213.306642558 |
| Directional search reduction | 0.957887335617 |
| Best parameter block | weights |
| Best block trial objective | 213.496613388 |
| Coefficient coordinate candidates | 64 |
| Best coefficient coordinate | absorption g25 negative step 1.0 |
| Best coefficient coordinate objective | 93.9257170886 |
| Greedy coefficient descent iterations | 28 |
| Greedy coefficient descent final objective | 8.84108251644 |
| Greedy coefficient descent source | cached checkpoint |
| Joint coefficient direction candidates | 12 |
| Best joint coefficient direction | nonzero_absorption_scales |
| Best joint coefficient objective | 8.82750258396 |
| Post-joint refinement iterations | 4 |
| Post-joint refinement final objective | 8.7132910678 |
| Post-coefficient weight refinement iterations | 1 |
| Post-coefficient weight refinement final objective | 8.7132910678 |
| Finite-difference coefficient-gradient candidates | 68 |
| Finite-difference coefficient-gradient norm | 31.9620610059 |
| Finite-difference coefficient-gradient best objective | 9.60321909946 |
| Finite-difference coefficient-gradient improved | false |
| Topology-neighbor candidates | 8 / 16 |
| Best topology-neighbor move | g9 -> g8 |
| Best topology-neighbor objective | 51.5218623108 |
| Topology scan status | all_32x16_topologies_fail_forcing_gate |
| Topology candidates | 15 |
| Best topology forcing lower bound | 8.61439354943 |
| Subset-search scan status | all_subset_topologies_fail_forcing_gate |
| Subset-search candidates | 4 |
| Best subset-search forcing lower bound | 23.296656694 |
| Hard threshold objective target | 1.0 |
| Final objective / target | 8.71329055076 |
| Final worst case | ecckd_rcemip_style_column_subset |
| Final worst metric | toa_forcing_max_abs |
| Final worst metric value | 2.61398716523 |
| Final worst metric threshold | 0.3 |
| Targeted worst-metric candidates | 96 |
| Targeted worst-metric objective | 8.71859705851 |
| Targeted full-objective candidate | 8.7132910678 |
| Targeted move accepted | false |
| Separated-component parameter count | 64 |
| Separated-component candidates | 64 |
| Separated-component target objective | 8.72141867832 |
| Separated-component full-objective candidate | 8.7132910678 |
| Separated-component move accepted | false |
| Table-refinement iterations | 2 |
| Table-refinement accepted moves | 2 |
| Table-refinement final objective | 8.71329055076 |
| Table-refinement improved | true |
| Next optimizer target case | ecckd_rcemip_style_column_subset |
| Next optimizer target metric | toa_forcing_max_abs |
| Required absolute reduction | 2.31398716523 |
| Required relative reduction | 0.885233 |
| Warm-start topology candidates | 2 / 16 |
| Warm-start topology best move | g4 -> g3 |
| Warm-start topology best objective | 8.7132910678 |
| Warm-start topology improved | false |
| Ranked topology candidates | 4 / 16 |
| Ranked topology best move | g21 -> g20 |
| Ranked topology best initial objective | 11.4112325601 |
| Ranked topology best refined objective | 8.7132910678 |
| Ranked topology improved | false |

The parameterization optimizes 16 shortwave spectral weights, 16 absorption coefficient scales, and 16 Rayleigh coefficient scales on the official ecCKD reduced hard-gate cases.

The topology scan ranks the existing official 32x16 reduced-accuracy candidates by the forcing-error lower bound implied by their TOA and surface forcing errors. This keeps the next optimizer tied to measured official-case evidence instead of only the current weighted-greedy topology.

Next required work: Move beyond bounded pressure-band table scales: run a stronger joint coefficient/table optimizer against flux and heating residuals or jointly optimize the reduced quadrature definition; the current table-refined 48-parameter path remains far above the hard-gate target.

Reactant status: `skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS`

Reactant surrogate check: `skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS`

Enzyme status: `skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS`

Enzyme surrogate check: `skipped_by_RH_SKIP_OPTIONAL_AD_CHECKS`
