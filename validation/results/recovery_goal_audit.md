# Recovery Goal Audit

Status: **not_complete**

- Blocked requirements: 0
- Partial requirements: 3
- Unmet requirements: 3
- Teacher-student recovery status: `passed`
- Objective reconstruction status: `ready_to_reconstruct_original_objective`
- Official training artifact status: `partial`
- Official training final objective / target: `8.605003990705882`
- Training recovery target status: `partial`
- Training target reduced scheme objective / target: `0.9994223340489347`
- CKDMIP preflight status: `ready_for_original_ecckd_objective`
- Derived flux generation plan status: `ready_for_original_ecckd_objective`
- Original objective assets ready: true
- Derived flux products: 18/18 final present, 0 with raw chunks
- Derived raw chunks: 0/90
- Completed-equivalent derived raw chunks: 90/90
- Observed derived raw chunk rate: `nothing` chunks/hour
- Estimated derived raw chunk hours remaining: `nothing`
- `ncrcat` concat tool: `present=true, path=/shared/home/greg/.local/bin/ncrcat, julia_concat_shim=true`
- Pareto points: 106
- Hard boundary forcing threshold: `0.3` W m^-2
- Official 32x32 worst boundary forcing error: `0.014033547037797689` W m^-2
- Published model accuracy status: `failed_threshold` (1/3 passing)
- Reduced near-miss limiter: `heating_rate_rmse` at `12.33245448795476`x threshold

## Requirements

| Requirement | Status | Finding | Evidence |
|---|---|---|---|
| Demonstrate parity with ecRad for full-accuracy models and reduced ecCKD models such as 16- and 32-band variants. | partial | The 32x32 published ecCKD model and at least one reduced candidate pass the hard thresholds; the published-model accuracy diagnostic currently records 1/3 published full-accuracy combinations passing the clean package-native reference gate. | `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecrad_accuracy_gate.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecrad_all_sky_ifs_gate.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_published_model_accuracy.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_accuracy.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_band_accuracy_pareto.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_refit_breakdown.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_scan.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_descent.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_descent_continuation.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_boundary_polish.json` |
| Demonstrate that reduced models meet accuracy criteria when compared to RRTMGP for representative atmosphere states. | partial | RRTMGP comparison metrics are emitted and at least one reduced candidate passes hard thresholds. | `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_32g_rrtmgp_comparison.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_accuracy.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_band_accuracy_pareto.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_refit_breakdown.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_scan.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_descent.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_descent_continuation.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/reduced_ecckd_leave_one_out_weight_coordinate_boundary_polish.json` |
| Demonstrate that the new radiation models can be integrated into Breeze simulations dynamically without blowing them up. | passed | The dedicated Breeze RCEMIP-style H100 artifact records a supported official ecCKD 32/32 runtime, passed gas-model accuracy status, and a finite >=4x RRTMGP speedup. | `/shared/home/greg/Projects/BreezeRadiativeHeatingDev/Breeze.jl/benchmarking/results/rcemip_h100_32x32x64/radiative_heating_rcemip_latest.json` |
| Reimplement the ecCKD training pipeline and demonstrate success by recovering one published model while varying only optimizer settings. | partial | Reactant/Enzyme checks pass and the current official/reduced optimizer reduces the objective, but the final objective remains 8.605003990705882x the hard target; CKDMIP assets are ready for original-objective recovery. | `/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_teacher_student_recovery_scan.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/official_ecckd_training.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_training_recovery_targets.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_objective_reconstruction_check.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ckdmip_training_data_preflight.json`<br>`/shared/home/greg/Projects/AnalyticBandRadiation.jl/validation/results/ecckd_derived_flux_generation_plan.json` |

## Prompt-to-Artifact Checklist

| Requirement ID | Covered | Status | Evidence Count | Gap |
|---|---:|---|---:|---|
| `ecrad_full_and_reduced_parity` | false | partial | 10 | The 32x32 published ecCKD model and at least one reduced candidate pass the hard thresholds; the published-model accuracy diagnostic currently records 1/3 published full-accuracy combinations passing the clean package-native reference gate. |
| `reduced_vs_rrtmgp_representative_states` | false | partial | 8 | RRTMGP comparison metrics are emitted and at least one reduced candidate passes hard thresholds. |
| `breeze_dynamic_integration` | true | passed | 1 |  |
| `reactant_enzyme_ecckd_training_recovery` | false | partial | 6 | Reactant/Enzyme checks pass and the current official/reduced optimizer reduces the objective, but the final objective remains 8.605003990705882x the hard target; CKDMIP assets are ready for original-objective recovery. |

## Quantitative Reduced-Model Status

- Best reduced candidate so far: 32x31 (63 total g-points), passed=true.
- Worst boundary forcing error: 0.2998267002146804 W m^-2 (TOA 0.2998267002146804, surface 0.2781916122328312).
- Method: official ecCKD 32x31 leave-one-out g23 support with exact boundary-polished quadrature weights
- Boundary-best dense leave-one-out diagnostic: 32x31, omitted SW g-point 25, status=failed_threshold.
- Near-miss limiter: ecckd_clear_sky_tropical_column heating_rate_rmse = 0.616622724397738 / 0.05 = 12.33245448795476x.
- Objective-best dense leave-one-out diagnostic: omitted SW g-point 23, objective=1.624423143786563x, limiter heating_rate_rmse = 0.08122115718932815 / 0.05 = 1.624423143786563x.
- Exact weight-coordinate improvement: omitted SW g-point 23, accepted=true, objective 1.624423143786563 -> 1.485564570060626, boundary 0.2547478587356409 W m^-2, heating RMSE 0.0742782285030313 K day^-1.
- Exact weight-coordinate descent: omitted SW g-point 23, accepted moves=6, objective 1.624423143786563 -> 1.2663455530591998, boundary 0.29987447424105085 W m^-2, heating RMSE 0.06331727765296 K day^-1.
- Exact weight-coordinate descent continuation: omitted SW g-point 23, accepted moves=8, objective 1.2663455530591998 -> 1.083385303619408, boundary 0.2994334374521941 W m^-2, heating RMSE 0.054169265180970406 K day^-1.
- Exact weight-coordinate boundary polish: omitted SW g-point 23, accepted moves=17, passed=true, objective 1.083385303619408 -> 0.9994223340489347, boundary 0.2998267002146804 W m^-2, heating RMSE 0.04987699106081957 K day^-1.

## Quantitative Training-Recovery Status

- Artifact status: partial; optimizer: deterministic multi-stage reduced ecCKD optimizer chain; parameters: 48; trainable SW g-points: 16.
- Objective: initial 214.26452989359106, final 8.605003990705882, target 1.0, final/target 8.605003990705882.
- Reactant check: passed; Enzyme check: passed; hard accuracy target met: false.
- Gap status: far_above_objective_target. Next work: Move beyond bounded pressure-band table scales: run a stronger joint coefficient/table optimizer against flux and heating residuals or jointly optimize the reduced quadrature definition; the current table-refined 48-parameter path remains far above the hard-gate target.

## Remaining Required Work

The CKDMIP data and derived ecCKD training flux products are ready for original-objective reconstruction. The remaining required work is running and validating the Reactant/Enzyme optimizer against the fixed published objective until one published model is recovered quantitatively.
