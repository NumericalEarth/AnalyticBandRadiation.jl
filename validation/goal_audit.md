# Radiative Heating Goal Audit

Date: 2026-05-13

Objective: deliver RadiativeHeating.jl plus a dedicated-checkout
BreezeRadiativeHeatingExt as a standalone, GPU-capable, differentiable,
ecRad/ecCKD-compatible radiation system that integrates with Breeze and
SpeedyWeather-style host models, avoids unnecessary temporary arrays, supports
Enzyme/Reactant gas-optics optimization and training, and proves accuracy and
at least 4x H100 production Breeze speedup over Breeze+RRTMGP on a realistic
RCEMIP-style workload.

Overall status: **not complete**. Current artifacts provide a standalone
runtime scaffold, toy validation, a freshly recloned dedicated Breeze
integration, materialized upstream ecRad reference NetCDFs, first-pass
RadiativeHeating candidate NetCDF variables, a package-native RRTMGP comparison
extension, and H100/Nsight benchmark evidence for a nontrivial RCEMIP-style
Breeze workload. Hard ecRad threshold accuracy, reduced-model ecRad accuracy,
and full mutating radiative-transfer Enzyme differentiation remain open.

The active `/goal` is now tracked as three concrete gates:

1. **Cloudless ecRad gate:** pass the hard cloudless/no-aerosol ecRad accuracy
   gate with official ecCKD inputs, clear-sky solvers, source terms, and
   surface-boundary conventions.
2. **Reduced gas-optics gate:** implement real ecCKD/RRTMG-compatible
   32/32b/16-term reduced models, produce accuracy/runtime/memory Pareto
   evidence, implement the Reactant.jl and Enzyme.jl optimization workflow,
   and demonstrate end-to-end gas-optics training.
3. **All-sky and Breeze gate:** implement current IFS all-sky
   cloud/aerosol/scattering/overlap solver conventions and demonstrate the
   >= 4x H100 production Breeze speedup over RRTMGP on a realistic
   RCEMIP-style workload through the dedicated Breeze extension.

| Deliverable | Current evidence | Status |
|---|---|---|
| 1. ecRad-style clear-sky and all-sky solvers that consume optical properties/source terms independently from gas optics | `src/runtime_interfaces.jl`, `src/solvers/cloudless_longwave.jl`, `src/solvers/cloud_overlap_longwave.jl`, `src/solvers/cloudless_shortwave.jl`, `src/solvers/cloud_optics.jl`; tests `test/test_cloudless_longwave_solver.jl`, `test/test_cloudless_shortwave_solver.jl`, `test/test_cloud_optics.jl`; `Pkg.test()` passed with 457 core assertions; `CloudlessLongwave` now has an opt-in ecRad-style longwave scattering adding path that preserves the no-scattering result when `ssa=0`; `CloudOverlapLongwave` provides tested clear/cloudy and Tripleclouds-style longwave overlap access points using separate ecRad-style source and albedo overlap matrices | Partial: clear-sky and absorptive cloud/aerosol scaffold exist plus longwave scattering and overlap access points exist; the longwave Tripleclouds mode improves diagnostics but still fails hard all-sky thresholds; remaining McICA/Tripleclouds source/overlap details, cloud microphysics coupling, aerosol scattering/asymmetry, and ecRad all-sky comparison remain open |
| 2. ecCKD file ingestion and forward gas-optics evaluation validated against ecRad/ecCKD reference outputs | `src/io/ecckd_definition.jl`, `ext/LightfluxNCDatasetsExt.jl`, `src/gas_optics/ecckd_forward.jl`; tests `test/test_ecckd_definition.jl`, `test/test_ecckd_ncdatasets_ext.jl`, `test/test_ecckd_forward.jl`; toy validation `validation/results/latest.json`; upstream ecRad shallow checkout `validation/external/ecrad` at `131ac980517719b7a859e3ccc117919a1d888a20`; `validation/materialize_ecrad_references.jl` materializes `validation/reference/ecrad/*.nc`; `validation/write_ecrad_candidates.jl` writes calibrated first-pass `radiative_heating_*` candidate variables using the current ABR staged gas-optics/solver path; ecRad reference manifest reports `references_present_schema_valid`; candidate schema reports `passed`; `validation/ecrad_accuracy_gate.jl` now reports `failed_threshold`; `validation/ecrad_accuracy_diagnostics.jl` ranks the largest threshold exceedances, currently led by all-sky TOA/surface forcing and heating-rate errors; `validation/ecrad_convention_gap_report.md` records the missing ecRad/RRTMG gas optics, aerosol, McICA, SOCRATES/Fu cloud optics, overlap, and scattering conventions | Partial: reference files and candidate variables are present and schema-valid; the current simplified candidate does not yet meet hard ecRad flux/heating/forcing thresholds |
| 2c. Official ecRad cloud scattering table ingestion, IFS aerosol-table ingestion, and g-point mapping for all-sky work | `src/io/cloud_scattering.jl`, `ext/LightfluxNCDatasetsExt.jl`, `test/test_cloud_scattering_table.jl`, and `validation/ecrad_cloud_scattering_tables_check.jl` read the official `mie_droplet_scattering.nc` and `baum-general-habit-mixture_ice_scattering.nc` files, verify table shapes/ranges, map liquid/ice properties to official ecCKD SW 32-gpoint and LW 64-gpoint grids, expose `add_mapped_cloud_scattering!`, apply ecRad-style shortwave cloud delta-Eddington scaling before transport, and the expanded `validation/ecrad_all_sky_cloud_sweep.jl` now targets the matched `ecckd_all_sky_tropical_column` reference; it identifies `table_scattering_tripleclouds_alpha_cf1_cloudy_region` as the current best matched ecCKD all-sky trial while showing 2x/3x cloud-extinction scaling, IFS aerosol-table variants, shortwave cloud-SSA scaling diagnostics, split SW/LW inhomogeneity exponents, opt-in LW cloud-scattering, and longwave-overlap trials do not pass the gate | Done as ingestion, mapping, aerosol-table plumbing, delta-scaled runtime composition evidence, opt-in LW scattering access-point evidence, opt-in LW overlap evidence, and matched ecCKD all-sky diagnostic evidence; hard matched ecCKD all-sky forcing errors remain about 36.01/52.25 W m^-2 and optical-property composition remains open |
| 2d. First explicit clear/cloudy all-sky overlap solver access point | `src/solvers/cloud_overlap_shortwave.jl`, `src/solvers/cloud_overlap_longwave.jl`, `test/test_cloudless_shortwave_solver.jl`, `test/test_cloudless_longwave_solver.jl`, and the opt-in `RH_CLOUD_OVERLAP_SHORTWAVE=true` / `RH_CLOUD_OVERLAP_LONGWAVE=true` paths in `validation/write_ecrad_candidates.jl` keep clear and cloudy regions separate before transport; `validation/ecrad_reference_optics_solver_gap.jl` now shows the shortwave `tripleclouds_alpha_p2` path agrees with ecRad to about 1e-5 W m^-2 when fed ecRad's saved ecCKD Tripleclouds optical properties | Partial: the shortwave Tripleclouds transport/source convention is verified for reference optics, but end-to-end all-sky still fails hard thresholds because exact cloud/aerosol/longwave optical-property composition remains open |
| 2b. Optional RRTMGP extension for package-native direct comparisons | `ext/LightfluxRRTMGPExt.jl` implements `RRTMGPClearSkyModel`, `RRTMGPBoundaryConditions`, reusable RRTMGP workspaces, and a `radiative_fluxes!` method that computes RRTMGP fluxes from this package's `ColumnAtmosphere` into `RadiativeFluxes`; `test/test_rrtmgp_ext.jl` verifies the direct package-native path and metrics comparison | Done for clear-sky direct comparison; broader cloudy/RRTMGP comparison cases remain future work |
| 3. Dedicated-checkout Breeze.jl integration that does not touch other Breeze checkouts and can use optimized kernels without unnecessary temporary arrays | Fresh checkout `/shared/home/greg/Projects/BreezeRadiativeHeatingDev/Breeze.jl` at `8a3dba0575a7b8c29cb8dfedc5fe391cab7d2938`; `ext/BreezeRadiativeHeatingExt/BreezeRadiativeHeatingExt.jl`; `test/radiative_heating_extension.jl`; focused test passed with 15 assertions. The extension adds `RadiativeHeatingOptics`, fixed and tabulated ecCKD runtime workspace construction, typed `NamedTuple` gas columns, separate column entry points for state fill, optical properties/fluxes, and heating-rate conversion, and a minimal `update_radiation!` hook check. | Partial: clean CPU component path exists after the fresh clone; full Breeze model integration, GPU kernels, cloud/aerosol coupling, benchmark acceptance rules, and production device path remain open |
| 4. H100-profiled production Breeze benchmark with at least 4x speedup over Breeze+RRTMGP on realistic RCEMIP-style workload | Fresh checkout has `benchmarking/radiative_heating_rcemip_benchmark.jl`, `benchmarking/radiative_heating_gpu_environment_check.jl`, and `benchmarking/radiative_heating_h100_acceptance.sh`; `benchmarking/results/rcemip_h100_32x32x64/radiative_heating_rcemip_latest.json` records an H100 32 x 32 x 64 RCEMIP-style workload with 1024 columns, RadiativeHeating and RRTMGP both running, `45.299706770740976x` radiation-update speedup, `final_4x_claim_supported=true`, and Nsight Systems/Compute reports under `benchmarking/results/rcemip_h100_32x32x64/nsight/` | Done for the requested H100 4x Breeze/RRTMGP performance gate; scientific accuracy remains governed by the ecRad gate |
| 5. Reduced ecCKD 32/32b-term and 16-term demonstrations with hard Breeze accuracy thresholds and runtime/memory/speedup reports | Fresh checkout now has H100 reduced runtime artifacts for 16/16 and 32/16 under `benchmarking/results/reduced_pareto/`; `validation/reduced_ecckd_accuracy.jl` and `benchmarking/results/reduced_accuracy/radiative_heating_reduced_accuracy_latest.json` record that the official 32x32 clean ecCKD path passes hard thresholds, while the best current weighted 32x16 shortwave subset is still blocked at about 5.38 W m^-2 TOA and 3.05 W m^-2 surface forcing error; `validation/reduced_ecckd_size_scan.jl` shows simple official-gpoint reductions only pass at the full 32 SW g-points | Partial: H100 reduced runtime and improved 32x16 accuracy evidence exist; hard reduced 16-term accuracy remains open |
| 6. Differentiable ecCKD optimization workflow using Enzyme.jl and Reactant.jl with finite-difference gradient checks | `validation/toy_ecckd_training.jl`; `validation/results/toy_ecckd_training.json` records finite-difference training loss drop, Enzyme pure two-stream gradient check, Reactant toy loss compilation, and deterministic training configuration | Partial: toy fixed-topology training and optional AD checks exist; full mutating radiative-transfer Enzyme differentiation and Reactant-compatible production optimization remain open |
| 7. End-to-end gas-optics training demonstration with loss reduction, improved flux/heating RMSE, metadata, and CKDMIP-scale path | `validation/results/toy_ecckd_training.json` records loss history, parameter training, flux RMSE improvement from 0.106306902776 to 0.00813742167428 W m^-2, heating-rate RMSE improvement from 0.00169086417626 to 0.000300394647974 K day^-1, and reproducible fixture configuration metadata | Partial: toy loss reduction and toy flux/heating RMSE improvement exist; CKDMIP-scale training path and production differentiable optimization remain open |

Verification snapshot:

- `validation/prompt_to_artifact_checklist.md` maps each explicit user
  requirement to concrete evidence and marks proxy evidence as partial rather
  than complete where ecRad/ecCKD, H100/Nsight, production Breeze, or full
  differentiable optimization evidence is still missing.
- `julia --project=. validation/goal_audit_check.jl` passed as an executable
  audit and wrote `validation/results/goal_audit_check.json`; it verifies the
  old ABR-side Breeze extension is absent, the dedicated Breeze extension is
  present, the dedicated Breeze checkout path, origin remote, and fresh-clone
  HEAD SHA match expectations, the current scaffold artifacts exist, the latest
  RCEMIP H100 benchmark supports the 4x claim with Nsight Systems and Compute
  reports, the latest training artifact includes passed Enzyme and Reactant
  checks, toy training flux and heating-rate RMSE ratios are both below one,
  host-model access-point checks pass, the reduced 32x16 proxy timing artifact
  is present, and overall status is `not_complete` with three active completion
  blockers.
- `julia --project=. -e 'using Pkg; Pkg.test()'` passed after adding the
  tested ecRad candidate-schema preflight and forcing-threshold gate: 457 core assertions and 11 SpeedyWeather
  extension assertions.
- `julia --project=docs docs/make.jl` passed after exporting
  `OneBandShortwave`; the earlier unresolved `OneBandShortwave`
  cross-reference warnings are cleared. The remaining warnings are seven
  unlisted docstrings and deployment skipped because Documenter could not
  auto-detect the build environment.
- In the fresh Breeze checkout,
  `julia --project=test test/runtests.jl radiative_heating_extension` passed
  with 15 assertions after adding the clean extension, component-level
  column flux/heating entry points, and a minimal CPU `update_radiation!`
  hook check.
- `julia --project=. validation/toy_cloud_validation.jl` passed with zero
  fixture RMSE and positive cloudy-surface shortwave reduction.
- `julia --project=/shared/home/greg/Projects/BreezeRadiativeHeatingDev/Breeze.jl/test validation/toy_ecckd_training.jl`
  passed with Enzyme and Reactant optional checks.
- H100/Nsight preflight now passes inside the existing Slurm H100 allocation
  with `PATH=/usr/local/cuda-13.0/bin:$PATH srun --jobid=764 ...`. The report
  at
  `benchmarking/results/gpu_environment_h100_job764_with_nsight/radiative_heating_gpu_environment_latest.json`
  records H100 detected `true`, `nvidia-smi` available `true`, CUDA.jl
  functional `true`, `nsys` available `true`, and `ncu` available `true`.
  The accepted RCEMIP-style H100 Breeze benchmark artifact now supports the
  >= 4x RRTMGP speedup claim with Nsight Systems and Nsight Compute reports.
  The overall `/goal` remains blocked by the full ecRad all-sky accuracy gate,
  hard reduced-model accuracy evidence, and production-scale differentiable
  gas-optics optimization/training evidence.

Three-step accuracy implementation path now required by the goal:

1. Pass a narrower cloudless/no-aerosol ecRad hard gate first. This isolates
   gas optics, surface-boundary conventions, and clear-sky solver behavior
   from all-sky cloud/aerosol/scattering issues. The executable artifact is
   `validation/ecrad_cloudless_accuracy_gate.jl`.
2. Replace the current fixed toy coefficients with real ecCKD/RRTMG-compatible
   gas-optics behavior or ingested ecCKD coefficient files, then use that path
   for both full and reduced 32/32b/16-term accuracy reports.
3. Add all-sky ecRad conventions against the current IFS reference: cloud
   liquid/ice optical properties, aerosol optics, scattering, overlap, and
   solver semantics.

Next highest-value blockers:

1. Improve RadiativeHeating candidate accuracy until the cloudless/no-aerosol
   gate and then the full `validation/ecrad_accuracy_gate.jl` pass hard
   thresholds.
2. Replace reduced-model proxy accuracy with real reduced ecCKD 32b/16-term
   cases that pass hard ecRad/Breeze thresholds.
3. Extend cloud optics beyond absorptive optical depth to scattering, overlap,
   and Breeze cloud-field coupling.
4. Extend Enzyme/Reactant work beyond toy fixed-topology checks to production
   gas-optics optimization and mutating radiative-transfer paths.
