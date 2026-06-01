# ecCKD Original Objective Terms

Status: **objective_terms_captured**

These terms and pass settings are fixed inputs for the Reactant/Enzyme original-objective recovery; optimizer settings may vary, but source data, loss terms, gases, weights, band mappings, and pass order may not.

Implementation status: `terms_captured_not_yet_recovered`

## Fixed Pass Sequences

| Kind | Application | Band structure | Tolerances | Common options | Modes |
|---|---|---|---|---|---|
| longwave | climate | fsck | 0.061, 0.0161 | `COMMON_OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_profile_weight=0.2 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 spectral_boundary_weight=0.1"` | relative-base -> relative-ch4 -> relative-n2o -> relative-cfc |
| shortwave | climate | rgb | 0.16, 0.047 | `COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.4 flux_profile_weight=0.1 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000 bounded_optimization=0"` | relative-base -> relative-ch4 -> relative-n2o |

## Objective Terms

| Kind | Term | Present | Source line | Description |
|---|---|---:|---:|---|
| longwave | `spectral_band_heating_rate` | true | 196 | Per-band squared heating-rate residual weighted by layer_weight and 86400^2. |
| longwave | `spectral_surface_down_flux` | true | 198 | Per-band surface downwelling flux residual weighted by flux_weight. |
| longwave | `spectral_toa_up_flux` | true | 199 | Per-band TOA upwelling flux residual weighted by flux_weight. |
| longwave | `interior_flux_profile` | true | 201 | Interior upwelling and downwelling profile residuals weighted by flux_profile_weight and interface weights. |
| longwave | `broadband_heating_rate` | true | 210 | Broadband heating-rate residual formed after summing over bands and blended by broadband_weight. |
| longwave | `broadband_boundary_flux` | true | 212 | Broadband surface-down and TOA-up residuals formed after summing over bands. |
| longwave | `spectral_boundary_flux` | true | 129 | Optional g-point surface-down and TOA-up boundary residuals weighted by spectral_boundary_weight. |
| longwave | `relative_reference_flux_subtraction` | true | 163 | Optional subtraction of relative-to CKD fluxes before evaluating residuals. |
| shortwave | `shortwave_transfer_branch` | true | 157 | Direct shortwave transfer for zero albedo; no-Rayleigh two-stream branch for reflective surfaces. |
| shortwave | `downwelling_only_heating_rate` | true | 197 | Heating-rate residual is computed from downwelling fluxes only in the CKD objective. |
| shortwave | `spectral_band_heating_rate` | true | 211 | Per-band squared heating-rate residual weighted by layer_weight and 86400^2. |
| shortwave | `spectral_surface_down_flux` | true | 213 | Per-band surface downwelling flux residual weighted by flux_weight. |
| shortwave | `spectral_toa_up_flux_20x` | true | 214 | Per-band TOA upwelling residual weighted by flux_weight with the official 20x multiplier. |
| shortwave | `interior_flux_profile` | true | 221 | Interior upwelling and downwelling profile residuals weighted by flux_profile_weight and interface weights. |
| shortwave | `broadband_heating_rate` | true | 247 | Broadband heating-rate residual formed after summing over bands and blended by broadband_weight. |
| shortwave | `broadband_surface_down_flux` | true | 250 | Broadband surface-down residual formed after summing over bands. |
| shortwave | `spectral_boundary_surface_down_flux` | true | 273 | Optional g-point surface-down residual weighted by spectral_boundary_weight. |
| shortwave | `relative_reference_flux_subtraction` | true | 164 | Optional subtraction of relative-to CKD fluxes before evaluating residuals. |

## Pass Details

### longwave

| Mode | Present | Gases | Training files | Specific options | Extra args |
|---|---:|---|---|---|---|
| `relative-base` | true | `composite h2o o3 co2` | `ckdmip_evaluation1_lw_fluxes_rel-180.h5  ckdmip_evaluation1_lw_fluxes_rel-280.h5 ckdmip_evaluation1_lw_fluxes_rel-415.h5  ckdmip_evaluation1_lw_fluxes_rel-560.h5 ckdmip_evaluation1_lw_fluxes_rel-1120.h5 ckdmip_evaluation1_lw_fluxes_rel-2240.h5` | `nothing` | `nothing` |
| `relative-ch4` | true | `ch4` | `ckdmip_evaluation1_lw_fluxes_present.h5  ckdmip_evaluation1_lw_fluxes_ch4-350.h5 ckdmip_evaluation1_lw_fluxes_ch4-700.h5  ckdmip_evaluation1_lw_fluxes_ch4-1200.h5 ckdmip_evaluation1_lw_fluxes_ch4-2600.h5 ckdmip_evaluation1_lw_fluxes_ch4-3500.h5` | `convergence_criterion=0.0005 flux_weight=0.5 negative_od_penalty=1.0e1` | `$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5` |
| `relative-n2o` | true | `n2o` | `ckdmip_evaluation1_lw_fluxes_present.h5 ckdmip_evaluation1_lw_fluxes_n2o-190.h5    ckdmip_evaluation1_lw_fluxes_n2o-270.h5 ckdmip_evaluation1_lw_fluxes_n2o-405.h5    ckdmip_evaluation1_lw_fluxes_n2o-540.h5` | `convergence_criterion=0.0005 flux_weight=0.5` | `$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5` |
| `relative-cfc` | true | `cfc11 cfc12` | `ckdmip_evaluation1_lw_fluxes_present.h5 ckdmip_evaluation1_lw_fluxes_cfc11-0.h5    ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5 ckdmip_evaluation1_lw_fluxes_cfc12-0.h5    ckdmip_evaluation1_lw_fluxes_cfc12-550.h5` | `convergence_criterion=0.0005 flux_weight=0.2 remove_min_max=1` | `$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5` |

### shortwave

| Mode | Present | Gases | Training files | Specific options | Extra args |
|---|---:|---|---|---|---|
| `relative-base` | true | `composite h2o o3 co2` | `ckdmip_evaluation1_sw_fluxes_rel-180.h5  ckdmip_evaluation1_sw_fluxes_rel-280.h5 ckdmip_evaluation1_sw_fluxes_rel-415.h5  ckdmip_evaluation1_sw_fluxes_rel-560.h5 ckdmip_evaluation1_sw_fluxes_rel-1120.h5 ckdmip_evaluation1_sw_fluxes_rel-2240.h5` | `convergence_criterion=0.01 spectral_boundary_weight=0.0` | `nothing` |
| `relative-ch4` | true | `ch4` | `ckdmip_evaluation1_sw_fluxes_present.h5  ckdmip_evaluation1_sw_fluxes_ch4-350.h5 ckdmip_evaluation1_sw_fluxes_ch4-700.h5  ckdmip_evaluation1_sw_fluxes_ch4-1200.h5 ckdmip_evaluation1_sw_fluxes_ch4-2600.h5 ckdmip_evaluation1_sw_fluxes_ch4-3500.h5` | `convergence_criterion=0.0005 max_no_rayleigh_wavenumber=15000` | `$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_rel-415.h5` |
| `relative-n2o` | true | `n2o` | `ckdmip_evaluation1_sw_fluxes_present.h5 ckdmip_evaluation1_sw_fluxes_n2o-190.h5    ckdmip_evaluation1_sw_fluxes_n2o-270.h5 ckdmip_evaluation1_sw_fluxes_n2o-405.h5    ckdmip_evaluation1_sw_fluxes_n2o-540.h5` | `convergence_criterion=0.0005 max_no_rayleigh_wavenumber=15000 remove_min_max=1` | `$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_rel-415.h5` |

Next required work:
- Implement Julia/Reactant evaluation of calc_cost_function_ckd_lw and calc_cost_function_ckd_sw using these captured terms.
- Run optimizer-only recovery for one published target and compare recovered weights/tables against the published CKD definition.
- Promote the result into official_ecckd_training only after the quantitative recovery target is met.
