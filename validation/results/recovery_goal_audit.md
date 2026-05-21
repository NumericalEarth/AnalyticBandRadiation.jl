# Recovery Goal Audit

Status: **not_complete**

- Blocked requirements: 1
- Partial requirements: 2
- Unmet requirements: 3
- Teacher-student recovery status: `passed`
- Objective reconstruction status: `blocked_missing_original_training_assets`
- CKDMIP preflight status: `ready_for_derived_flux_generation`
- Derived flux generation plan status: `derived_flux_generation_required`
- Original objective assets ready: false
- Derived flux products: 3/18 final present, 1 with raw chunks
- Derived raw chunks: 3/90
- Completed-equivalent derived raw chunks: 18/90
- Observed derived raw chunk rate: `8.046269021356558` chunks/hour
- Estimated derived raw chunk hours remaining: `10.81246472981241`
- `ncrcat` concat tool: `present=true, path=/shared/home/greg/.local/bin/ncrcat, julia_concat_shim=true`
- Pareto points: 22
- Hard boundary forcing threshold: `0.3` W m^-2
- Official 32x32 worst boundary forcing error: `0.014033547037797689` W m^-2

## Requirements

| Requirement | Status | Finding | Evidence |
|---|---|---|---|
| Demonstrate parity with ecRad for full-accuracy models and reduced ecCKD models such as 16- and 32-band variants. | partial | Full official ecCKD 32x32 passes, but currently measured reduced 16-g shortwave candidates fail the hard thresholds. | `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecrad_accuracy_gate.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecrad_all_sky_ifs_gate.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_accuracy.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_band_accuracy_pareto.json` |
| Demonstrate that reduced models meet accuracy criteria when compared to RRTMGP for representative atmosphere states. | partial | RRTMGP comparison metrics are emitted for the official 32x32 path on representative states, but reduced 16-g candidates do not yet pass hard accuracy criteria. | `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_32g_rrtmgp_comparison.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_accuracy.json` |
| Demonstrate that the new radiation models can be integrated into Breeze simulations dynamically without blowing them up. | passed | The dedicated Breeze RCEMIP-style H100 artifact records a supported official ecCKD 32/32 runtime, passed gas-model accuracy status, and a finite >=4x RRTMGP speedup. | `/shared/home/greg/Projects/BreezeRadiativeHeatingDev/Breeze.jl/benchmarking/results/rcemip_h100_32x32x64/radiative_heating_rcemip_latest.json` |
| Reimplement the ecCKD training pipeline and demonstrate success by recovering one published model while varying only optimizer settings. | blocked | Teacher-student coefficient recovery passes and upstream CKDMIP data is present, but exact original-objective recovery is blocked until the derived ecCKD 5gas/rel training flux products are generated. | `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_teacher_student_recovery_scan.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_objective_reconstruction_check.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ckdmip_training_data_preflight.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_derived_flux_generation_plan.json` |

## Prompt-to-Artifact Checklist

| Requirement ID | Covered | Status | Evidence Count | Gap |
|---|---:|---|---:|---|
| `ecrad_full_and_reduced_parity` | false | partial | 4 | Full official ecCKD 32x32 passes, but currently measured reduced 16-g shortwave candidates fail the hard thresholds. |
| `reduced_vs_rrtmgp_representative_states` | false | partial | 2 | RRTMGP comparison metrics are emitted for the official 32x32 path on representative states, but reduced 16-g candidates do not yet pass hard accuracy criteria. |
| `breeze_dynamic_integration` | true | passed | 1 |  |
| `reactant_enzyme_ecckd_training_recovery` | false | blocked | 4 | Teacher-student coefficient recovery passes and upstream CKDMIP data is present, but exact original-objective recovery is blocked until the derived ecCKD 5gas/rel training flux products are generated. |

## Quantitative Reduced-Model Status

- Best reduced candidate so far: 32x16 (48 total g-points), passed=false.
- Worst boundary forcing error: 2.075241671208687 W m^-2 (TOA 2.075241671208687, surface 2.026956688318073).
- Method: weighted greedy 16 shortwave g-point subset with boundary-aware table, component, structural, objective-probe, surface-probe, capped table, continuation, post-capped weight, post-weight surface-table, bounded weight, four current component-scale refits, selected current gas-pressure component scan refit, gas-pressure continuation refit, weighted gas-pressure continuation refit, and high-weight gas-pressure continuation refit

## Current Blocker

Exact original ecCKD objective recovery now has the upstream CKDMIP tree available, but still requires local generation of the derived ecCKD `5gas-*` and `rel-*` training flux products. The current teacher-student recovery is useful AD evidence, but it is not a substitute for reconstructing the published objective.
