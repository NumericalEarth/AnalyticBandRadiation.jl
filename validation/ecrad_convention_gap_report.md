# ecRad Convention Gap Report

Date: 2026-05-13

This report records why the current `radiative_heating_*` candidate variables
do not pass `validation/ecrad_accuracy_gate.jl`. It is diagnostic only; the
hard pass/fail source of truth remains the accuracy gate.

## Reference Configuration

The materialized reference NetCDF files come from the upstream ecRad IFS test
data in `validation/external/ecrad/test/ifs`.

The all-sky reference log
`validation/external/ecrad/test/ifs/ecrad_meridian_default_out_REFERENCE.log`
records:

- Shortwave gas model: `RRTMG-IFS`
- Longwave gas model: `RRTMG-IFS`
- Aerosols: enabled, with general aerosol optics
- Clouds: enabled
- Liquid cloud optics: `SOCRATES`
- Ice cloud optics: `Fu-IFS`
- Cloud overlap: `Exp-Ran`
- Cloud PDF shape: `Gamma`
- Shortwave solver: `McICA`
- Longwave solver: `McICA`
- Longwave cloud scattering: enabled
- Cloud/aerosol/surface per-g-point properties: disabled for this RRTMG-IFS reference

The clear-sky reference used by
`validation/materialize_ecrad_references.jl` comes from
`ecrad_meridian_cloudless_out_REFERENCE.nc`, but still uses the ecRad/RRTMG
gas-optics and cloudless solver conventions rather than the simplified
fixed-coefficient candidate model.

## Current Candidate

`validation/write_ecrad_candidates.jl` writes candidate variables using:

- `EcCKDTabulatedGasOpticsModel` when `RH_CANDIDATE_GAS_OPTICS=official_ecckd`,
  using official ecCKD LW/SW coefficient files for H2O, O3, CO2, CH4, N2O,
  CFC11, and CFC12
- `CloudlessLongwave`
- `CloudlessShortwave` with a first conservative Rayleigh scattering
  approximation
- Broadband scalar surface albedo and surface longwave boundary matching
- Optional cloud optical-depth hook, disabled by default, with a
  phase-separated liquid/ice/cloud-fraction access point and separate
  shortwave absorption/scattering/asymmetry channels for all-sky plumbing
- Staged layer aerosol optics with separate shortwave
  absorption/scattering/asymmetry channels

This candidate is useful for plumbing validation because it writes the expected
variables with the correct shapes, but it is not yet physically equivalent to
the ecRad reference setup.

## Measured Gaps

`validation/ecrad_accuracy_diagnostics.jl` currently reports 36 failed metrics.
The largest threshold exceedances are:

- all-sky TOA forcing absolute error
- all-sky surface forcing absolute error
- all-sky shortwave flux errors
- all-sky and clear-sky heating-rate errors

The all-sky reference includes nonzero `cloud_fraction`,
`liquid_water_path`, `ice_water_path`, `re_liquid`, `re_ice`,
`inv_cloud_effective_size`, `fractional_std`, `overlap_param`, and
`aerosol_mmr`. The materializer now preserves the cloud effective-size,
inhomogeneity, overlap-parameter, and aerosol fields from the upstream ecRad
input so candidate physics can use the same information as the reference. The
candidate writer can pass cloud fraction and liquid/ice paths through a
phase-separated cloud-optics model and can partition shortwave extinction into
absorption and scattering with optical-depth-weighted asymmetry. The current
model does not yet use the preserved `fractional_std` and `overlap_param` in a
full Tripleclouds-equivalent way, and it is not yet an IFS-equivalent
SOCRATES/Fu-IFS scattering, phase-function, aerosol, and overlap treatment.

## Negative Result

A narrow test that added simple absorptive cloud optical depth from
`cloud_fraction * (liquid_water_path + ice_water_path)` made the all-sky
diagnostics worse for the trial coefficients tested. That result is consistent
with the ecRad reference requiring scattering, overlap, and solver convention
support rather than just extra absorption.

A later phase-separated scattering trial used explicit liquid/ice shortwave
extinction, high single-scattering albedo, and forward asymmetry. It reduced
the all-sky worst boundary forcing miss from about `572` W m^-2 to about
`217` W m^-2, but still failed the hard thresholds and worsened heating-rate
errors. This confirms that the new scattering channel is useful plumbing, but
that coefficient tuning alone is not a substitute for IFS-equivalent overlap,
cloud-source, and solver semantics.

A small cloud-fraction exponent sweep using the same scattering trial found
that `RH_CLOUD_FRACTION_EXPONENT=0.5` reduced the all-sky surface forcing miss
to about `165` W m^-2, but left the TOA forcing miss around `213` W m^-2. The
candidate writer records this exponent in
`radiative_heating_cloud_fraction_exponent`. This is another negative result:
deterministic cloud-fraction scaling is not enough to replace ecRad's overlap
and McICA-style all-sky semantics.

An opt-in effective-radius cloud-extinction path is available through
`RH_CLOUD_EFFECTIVE_RADIUS_OPTICS=true`. It uses the materialized `re_liquid`
and `re_ice` fields to estimate shortwave extinction as
`3 * water_path / (2 * density * effective_radius)`, with separate liquid and
ice scaling knobs. This path is still available as a diagnostic fallback, but
the best current all-sky sweep now uses official table-mapped cloud scattering
instead of the effective-radius extinction proxy.

An opt-in aerosol path is available through `RH_AEROSOL_OPTICS=true`. It
converts the materialized `aerosol_mmr` into a layer aerosol path and composes
it through the staged aerosol optical-depth API. Positive shortwave aerosol
mass-extinction trials on the current all-sky case worsened the boundary
forcing errors, so the current candidate artifact leaves aerosol optics off.
The path is retained as plumbing for future IFS aerosol optical-property work.

An official ecCKD ingestion trial is now available through
`RH_CANDIDATE_GAS_OPTICS=official_ecckd`. This mode reads the ecRad checkout's
LW 64-gpoint and SW 32-gpoint ecCKD definition files, preserves the H2O
mole-fraction axis for layer-dependent interpolation, converts H2O/O3 mass
mixing ratios and CO2/CH4/N2O volume mixing ratios to molar column amounts, and
writes normal `radiative_heating_*` candidate variables. Adding the official
longwave Planck source table improved the cloudless heating-rate RMSE from
about 69 K day^-1 to about 9.8 K day^-1.
Preserving the official pressure-dependent temperature grid, log-pressure
interpolation, ecRad's pressure-weighted full-level temperature convention, and
solar-path scaling from the materialized `cos_solar_zenith_angle` reduced some
cloudless misses. Dynamic H2O interpolation, Rayleigh scattering, the ecCKD
composite background-gas term, spectral surface emission, ecRad-style
shortwave direct/diffuse adding, and the ecRad no-scattering interface-source
formula, missing-coefficient zero filling for SW-only gaps, CFC11/CFC12
longwave absorption, and relative-linear reference mole fractions for CH4/N2O
closed the clean ecCKD cloudless/no-aerosol hard gate. The latest clean ecCKD
tropical TOA/surface mean biases are about `-0.028` and `-0.010` W m^-2; the
clean ecCKD RCEMIP-style subset is about `-0.027` and `-0.0045` W m^-2.

An opt-in longwave interface-source path is available through
`RH_LW_INTERFACE_SOURCES=true`. It uses ecRad's no-scattering linear Planck
source formula with half-level Planck values. It is now the default for
`RH_CANDIDATE_GAS_OPTICS=official_ecckd`, where the candidate writer also uses
the ecCKD Planck table to distribute the prescribed broadband surface longwave
flux spectrally. This combination reduced the clear-sky tropical longwave RMSEs
to about `3.34` W m^-2 for `lw_up` and `1.98` W m^-2 for `lw_down`.

The official ecCKD candidate now uses the gases present in the official
definition files and the ecRad IFS input: H2O, O3, CO2, CH4, N2O, CFC11, and
CFC12, plus the ecCKD composite background-gas absorption term. This is more
physically faithful than the earlier H2O/CO2-only candidate. The reader no
longer collapses H2O to one fixed mole-fraction slice, and it now applies the
official relative-linear reference mole fractions where present.

`validation/ecrad_flux_bias_diagnostics.jl` now writes signed flux-bias
diagnostics alongside the hard gate. The latest clean ecCKD candidate passes
the focused cloudless/no-aerosol gate. The full all-sky hard gate still fails,
so the next physics work is cloud/aerosol/scattering/overlap solver
conventions and reduced-model accuracy evidence.

`validation/ecrad_all_sky_cloud_effect_diagnostics.jl` now writes cloud-effect
diagnostics using the materialized ecRad clear-sky fluxes and matching
`radiative_heating_*_clear` candidate variables.
`validation/ecrad_all_sky_cloud_sweep.jl` now compares the older
effective-radius proxy, table-mapped SOCRATES/Fu cloud scattering, first
explicit overlap-blend trials, and the IFS aerosol-table path against the
matched `ecckd_all_sky_tropical_column` reference. This matters because the
older `all_sky_tropical_column` file is RRTMG-IFS + McICA, while the all-sky
candidate and saved radiative-property diagnostics are ecCKD + Tripleclouds.
The sweep leaves candidate NetCDF variables in the best tested state for the
matched ecCKD all-sky target. The current best trial is
`table_scattering_tripleclouds_alpha_cf1_cloudy_region`, with official ecCKD
gas optics, table-mapped cloudy-region shortwave cloud scattering,
ecRad-style shortwave cloud delta-Eddington scaling before transport, actual
cloud fractions in Tripleclouds transport, and the ecRad-configured
inhomogeneity-overlap exponent `2.0`. Its matched ecCKD all-sky TOA and
surface forcing misses are about `36.01` and `52.25` W m^-2, with worst
threshold ratio about `174.16`. This remains far outside the `0.3` W m^-2
hard boundary thresholds. Scaling the table-mapped liquid/ice cloud
extinction by factors of `2` and `3`, the current aerosol-table variants,
shortwave cloud-SSA scaling, split shortwave/longwave inhomogeneity exponents,
and opt-in longwave overlap/scattering trials all fail the gate. The matched
reference result localizes the remaining full-gate failure primarily to
all-sky optical-property generation/composition, especially exact aerosol and
longwave cloud conventions, not to the clean ecCKD cloudless gas-optics path
or the shortwave Tripleclouds transport when ecRad optical properties are
provided directly.

`validation/ecrad_all_sky_optics_gap.jl` now compares the candidate optical
properties to ecRad's saved all-sky `radiative_properties_ecckd_tc.nc` before
radiative-transfer solving using the current table-scattering best
configuration. The grid-mean candidate still has too little cloud optical
depth: `od_sw_cloud` mean is about `0.026` versus ecRad `0.256`, and
`od_lw_cloud` mean is about `0.031` versus ecRad `0.512`, because the current
runtime composition is still a grid-mean no-longwave-scattering shortcut. The
diagnostic now also maps the official liquid/ice cloud-scattering tables to the
32-gpoint LW and SW ecCKD grids with an ecRad-style triangular spectral mapper,
solar interval weights where available, delta-Eddington averaging, layer
effective radii, in-cloud condensate paths, and ecRad-style cloud
delta-Eddington scaling before transport. That cloudy-region view now matches
ecRad's saved shortwave cloud optical depth and asymmetry closely:
`od_sw_cloud` mean is about `0.255846` versus ecRad `0.255843`, and
`asymmetry_sw_cloud` mean is about `0.0950398` versus ecRad `0.0950405`.
The longwave cloudy-region optical depth is still high (`0.752` versus
`0.512`). This confirms the remaining shortwave mismatch is no longer a scalar
cloud-extinction calibration or nearest-midpoint mapping artifact: ecRad's
all-sky Tripleclouds/McICA path stores cloudy-region optical properties and
then uses cloud fraction, overlap, inhomogeneity, and solver semantics in the
transport. That is now the main all-sky convention gap.

The runtime API now exposes `CloudyRegionCloudOpticalProperties` and
`cloudy_region_optical_properties!` so host models and future all-sky solvers
can carry cloud fraction and overlap separately from cloudy-region optical
depth. This is only an interface/representation step: the optics-gap
diagnostic shows that simply removing cloud-fraction scaling does not match
ecRad by itself because the current provisional effective-radius optics still
uses simplified scattering properties (`ssa_sw_cloud=1` and large asymmetry)
rather than ecRad's spectrally averaged SOCRATES/Fu cloud optical tables.

The package now has a `CloudScatteringTable` API and a NCDatasets reader for
the official ecRad cloud-scattering inputs. The validation artifact
`validation/ecrad_cloud_scattering_tables_check.jl` verifies that the
`mie_droplet_scattering.nc` and
`baum-general-habit-mixture_ice_scattering.nc` files can be read, that their
coefficient arrays have the expected table shapes, and that sampled
mass-extinction, single-scattering-albedo, and asymmetry values are finite and
bounded. The same artifact now reads the ecCKD `gpoint_fraction` mappings and
projects both liquid and ice cloud-scattering tables to the official SW
32-gpoint and LW 64-gpoint grids with finite bounded properties. This closes
the table-ingestion and first g-point mapping parts of the all-sky optics work,
but not solver use.

An explicit `ShortwaveCloudOverlapOpticalProperties` /
`CloudOverlapShortwave` access point now keeps clear-region and cloudy-region
shortwave optical properties separate and blends interface fluxes with a
named cloud-fraction overlap rule. The ecRad candidate writer can exercise
this path with `RH_CLOUD_OVERLAP_SHORTWAVE=true`. The first deterministic
maximum-overlap trial is a negative result: with table-mapped cloud scattering
enabled, the all-sky tropical surface and TOA forcing misses increase to about
`637` and `651` W m^-2 in the expanded sweep. This confirms that simple clear/cloudy interface
blending is not an adequate substitute for ecRad's McICA/Tripleclouds
two-stream adding and cloud-region source semantics.

The candidate writer can also test the same overlap access point using true
cloudy-region shortwave optics, enabled by
`RH_CLOUD_OVERLAP_USE_CLOUDY_REGION_OPTICS=true`, where liquid and ice water
paths are converted from grid-mean to in-cloud paths before adding the official
table-mapped cloud scattering. This is also a negative result: the expanded
cloud sweep shows the cloudy-region average and maximum overlap trials still
produce TOA and surface forcing misses of about `619` and `634` W m^-2, with
heating cloud-effect errors around `814` to `1151` K day^-1. The failure is
therefore not just that the overlap blend used grid-mean cloud optical depth;
the missing mechanism is the ecRad adding method with overlap matrices,
inhomogeneous cloud-region fractions, and source coupling.

An intermediate `CloudOverlapShortwave(overlap = :adding)` mode now mixes
clear/cloudy layer reflectance and transmittance before the shortwave adding
pass rather than blending completed interface fluxes. In the all-sky sweep,
`table_scattering_adding_cloudy_region` reduces the extreme final-flux-blend
TOA/surface misses from about `619`/`634` W m^-2 to about `307`/`315` W m^-2,
but it is still much worse than the non-overlap `table_scattering` trial at
about `79`/`84` W m^-2. This confirms that adding during transport is necessary
but not sufficient; the remaining solver work needs region fractions,
inhomogeneity scaling, and overlap matrices analogous to ecRad's
Tripleclouds/McICA implementation.

A further `CloudOverlapShortwave(overlap = :matrix_maximum)` diagnostic carries
separate clear/cloudy region fluxes through a two-region maximum-overlap matrix
during the adding pass. This improves over scalar adding, reducing TOA/surface
misses to about `162`/`159` W m^-2 and the heating cloud-effect maximum to
about `10.9` K day^-1, but it still does not beat the current best
grid-mean `table_scattering` candidate. The ordering is now clear: final-flux
blending is worst, scalar adding is better, two-region matrix adding is better
again, and the remaining gap is likely the full ecRad three-region
Tripleclouds/McICA treatment with gamma-distribution inhomogeneity scaling,
overlap decorrelation, and cloud/source coupling.

The overlap container now also carries the materialized ecRad `overlap_param`,
and `CloudOverlapShortwave(overlap = :matrix_alpha)` uses it in a two-region
Hogan-Illingworth alpha-overlap matrix. In the current simplified two-region
solver this is a negative result: `table_scattering_matrix_alpha_cloudy_region`
misses TOA/surface forcing by about `183`/`184` W m^-2, worse than
`matrix_maximum` and still worse than the best grid-mean `table_scattering`
case. This localizes the next required implementation step to using
`fractional_std` to split the cloudy region into thin/thick Tripleclouds
regions with ecRad's gamma-distribution optical-depth scaling, rather than
only applying the overlap parameter to a single cloudy region.

`CloudOverlapShortwave(overlap = :tripleclouds_alpha)` now performs that first
three-region diagnostic split: each cloudy layer is divided into thin and
thick cloudy regions using ecRad's gamma optical-depth scaling from
`fractional_std`, and the cloudy-region transitions use an inhomogeneity
overlap derived from the supplied alpha overlap parameter. This improves over
the two-region alpha diagnostic but is still not the winning candidate:
`table_scattering_tripleclouds_alpha_cloudy_region` with the older
grid-mean-tuned cloud-fraction exponent `0.5` misses TOA/surface forcing by
about `161`/`156` W m^-2. This identified an important convention issue:
ecRad's overlap machinery uses the actual layer cloud fraction, not the
sqrt-cloud-fraction tuning that helped the provisional grid-mean shortcut.
When the Tripleclouds path uses `RH_CLOUD_FRACTION_EXPONENT=1.0`, the
TOA/surface forcing miss falls to about `32.23`/`50.75` W m^-2 at
inhomogeneity exponent `2.0`. The best tested variant uses exponent `4.0`,
with longwave Tripleclouds enabled and reaches about `46.86`/`47.35` W m^-2
and worst threshold ratio about `157.83`.
An explicit longwave cloud-scattering storage/solver access point now exists,
but the corresponding sweep trial worsens the current simplified transport
state (`39.99`/`52.19` W m^-2 TOA/surface), so it remains opt-in until it is
coupled to a full longwave Tripleclouds/McICA transport.
An explicit longwave clear/cloudy overlap access point also exists; it reduces
the heating cloud-effect maximum in the same sweep but does not improve the
hard boundary forcing gate (`42.31`/`50.06` W m^-2 TOA/surface), so it is also
kept out of the best candidate state.
The corrected longwave `tripleclouds_alpha` option now beats the two-region
longwave overlap diagnostic in the all-sky sweep, although it still falls far
outside the hard boundary forcing thresholds.
The remaining all-sky gap therefore is no longer cloud-region optical-depth
magnitude or cloud-fraction use; it likely needs the full ecRad
McICA/Tripleclouds source treatment, longwave three-region transport, aerosol
details, and exact cloud/source coupling.

`validation/ecrad_reference_optics_solver_gap.jl` now isolates shortwave solver
error by bypassing RadiativeHeating's gas/cloud/aerosol optical-property
generation and feeding the solvers ecRad's saved `radiative_properties_ecckd_tc`
shortwave optical properties, comparing against the matched
`ecckd_all_sky_tropical_column` fluxes. This diagnostic uses ecRad's
per-g-point solar weights and separate direct/diffuse surface albedos through
`ShortwaveBoundaryConditions(surface_albedo_direct = ...)`. With matched
reference optics, the clear-sky shortwave path is at roundoff-scale agreement
and `tripleclouds_alpha_p2` reaches about `1.4e-5`/`2.6e-5` W m^-2
TOA/surface net maximum error. This verifies the shortwave Tripleclouds
transport/source convention for the ecCKD Tripleclouds case when ecRad optical
properties are supplied directly. The remaining end-to-end all-sky failure is
therefore an optical-property generation/composition gap rather than a
shortwave Tripleclouds transport gap.

The shortwave two-stream direct-source clamp now distinguishes ecRad's two
helper conventions. Cloudless/homogeneous paths use the
`calc_reflectance_transmittance_sw` unit-normalized source bound, while
Tripleclouds-style paths use the `calc_ref_trans_sw` bound
`mu0 * (1 - trans_dir_dir)`. The reference-optics diagnostic is essentially
unchanged, so the remaining clear-sky solver gap is elsewhere in the
saved-flux convention or adding/source treatment.

The shortwave adding source expression has also been aligned with ecRad's
`adding_ica_sw`: the multiple-reflection denominator applies to the
transmitted diffuse/source term and not to the direct-reflection source term.
This is a small but real convention fix. In the end-to-end sweep it improved
the table-scattering-only all-sky TOA forcing miss to about `75.37` W m^-2,
while the surface forcing miss remained about `83.83` W m^-2.

The mapped cloud-scattering composition now also applies ecRad's shortwave
cloud delta-Eddington scaling before transport, matching the
`do_sw_delta_scaling_with_gases=false` IFS test configuration where cloud
optical depths are scaled before being combined with gas optical properties.
This reduced the best table-scattering all-sky forcing errors to about
`72.65` W m^-2 at TOA and `65.03` W m^-2 at the surface.

Using actual cloud fractions in the Tripleclouds overlap path, rather than the
grid-mean tuning exponent `0.5`, is another real convention fix. After adding
the matched ecCKD Tripleclouds all-sky reference, the reference-optics solver
gap shows that the ecRad-configured inhomogeneity-overlap exponent `2.0`
matches the saved ecRad shortwave fluxes when ecRad optical properties are
provided directly.

The materialized all-sky references now preserve ecRad's 32-gpoint diffuse and
direct shortwave surface albedos from `radiative_properties_ecckd_tc.nc`.
`ShortwaveBoundaryConditions` and the candidate writer can use separate
diffuse/direct albedos when the g-point data are available. This is a real
convention fix, but it is not sufficient by itself: after rerunning the cloud
sweep with cloud delta-Eddington scaling enabled and actual cloud fractions in
Tripleclouds transport against the matched ecCKD all-sky reference, the best
candidate reaches about `36.01`/`52.25` W m^-2 TOA/surface forcing error.

The latest optical-property diagnostic confirms that this is not mainly a
single scalar-tuning problem. Against ecRad's saved
`radiative_properties_ecckd_tc.nc`, the grid-mean candidate underestimates
cloud optical depth (`od_lw_cloud` mean about `0.030` versus ecRad `0.512`;
`od_sw_cloud` mean about `0.026` versus ecRad `0.256`), while the delta-scaled
cloudy-region candidate matches shortwave cloud optical depth
(`od_sw_cloud` mean about `0.255846`). A guarded
`RH_CLOUD_SCATTERING_SW_SSA_SCALE` diagnostic can reduce the remaining
cloud-SSA RMSE while preserving total extinction and asymmetry, but the
all-sky sweep shows `1.08` and `1.15` scaling worsen TOA forcing and do not
move the surface forcing enough to matter. The current blocker is therefore
not a scalar cloud-SSA tuning problem. A separate diagnostic split of shortwave
and longwave Tripleclouds inhomogeneity exponents also worsens the matched
ecCKD all-sky gate: `SW p8 + LW p4` and `SW p12 + LW p4` are both worse than
the ecRad-configured `p2` shortwave transport. The reference configuration also has
general aerosol optics enabled. The current best candidate leaves aerosol
optics off because both fixed-RH and diagnosed layer-RH IFS aerosol-table
variants still miss the all-sky hard thresholds by more than two orders of
magnitude and are slightly worse than the cloud-table-only trial after cloud
delta scaling. Closing the all-sky gate therefore requires an IFS-compatible
cloud-region optical-property path plus exact aerosol optical properties, not
only more end-to-end scale sweeps.

A first IFS aerosol-table candidate path has been added to the validation
writer. It maps `aerosol_ifs_49R1_20230119.nc` onto the ecCKD g-point grids,
uses the ecRad 49R1 aerosol type map, and composes per-type aerosol paths into
longwave absorption plus shortwave absorption/scattering/asymmetry. This is
now included in the all-sky cloud sweep as
`table_scattering_ifs_aerosol_table` and
`table_scattering_ifs_aerosol_table_fixed_rh08`. After adding cloud and
aerosol delta-Eddington scaling, the layer-RH lookup variant reaches about
`76.92` W m^-2 at TOA and `67.27` W m^-2 at the surface. The fixed-RH `0.8`
variant reaches about `76.41` W m^-2 TOA and `67.21` W m^-2 surface forcing
error. Both are slightly worse than the aerosol-free `table_scattering`
candidate, so the sweep leaves aerosol optics off for the current best state.
This is useful evidence that aerosols matter, but it is still not a passing
configuration.

The ecRad checkout contains `ecrad_meridian_ecckd_tc_out_REFERENCE.nc` and
`ecrad_meridian_ecckd_mcica_out_REFERENCE.nc`, including clear-sky flux
variables. Those are useful diagnostics, but they are not the hard
cloudless/no-aerosol target because they come from ecCKD all-sky configurations
with aerosol settings. A local `bin/ecrad` executable has now been built
against locally extracted NetCDF/HDF5 libraries, and the hard first gate now
includes clean `ECCKD + Cloudless + no aerosol` tropical and RCEMIP-style
reference files generated from
`test/ifs/ecrad_meridian_ecckd_cloudless_noaer_out.nc`.

## Required Implementation Work

Passing the hard ecRad gate requires real model work, not only coefficient
tuning. The `/goal` is now organized around three concrete gates, with the
all-sky work decomposed below:

1. Pass the cloudless/no-aerosol ecRad gate first. That gate isolates
   gas-optics, boundary, and clear-sky solver conventions from the all-sky
   physics and is written by `validation/ecrad_cloudless_accuracy_gate.jl`.
   This gate now passes for the clean ecCKD tropical and RCEMIP-style cases.
2. Implement ecRad/RRTMG-compatible gas-optics conventions or ingest real ecCKD
   coefficients for the reference cases, including the reduced 32/32b/16-term
   models needed for the Breeze Pareto reports.
3. Add current IFS all-sky conventions: liquid and ice cloud optical
   properties, effective radii or equivalent assumptions, aerosol optics,
   scattering, single-scattering albedo, asymmetry, cloudy-region optical
   depth storage, McICA/Tripleclouds-style overlap, and solver semantics for
   the selected reference case set. The
   runtime API now has a phase-separated liquid/ice/cloud-fraction access
   point and shortwave scattering/asymmetry channel, but it does not yet
   satisfy this gate. Aerosol optics now have the same staged
   absorption/scattering/asymmetry channels, and the official cloud
   scattering tables can now be ingested and mapped to the ecCKD g-point
   grids, and a first explicit clear/cloudy overlap access point exists, but
   the deterministic blend is not scientifically adequate. IFS-equivalent
   McICA/Tripleclouds transport, cloud-region source handling, inhomogeneity,
   and overlap semantics are still missing.
4. Re-run `validation/write_ecrad_candidates.jl`,
   `validation/ecrad_accuracy_gate.jl`, and
   `validation/ecrad_accuracy_diagnostics.jl` after each physics increment.
   For all-sky changes, also run
   `validation/ecrad_all_sky_cloud_effect_diagnostics.jl` and, when tuning
   provisional cloud optics, `validation/ecrad_all_sky_cloud_sweep.jl` to
   separate clear-sky drift from cloud radiative-effect errors.

Until those items are implemented, the goal remains blocked by the full all-sky
ecRad accuracy gate and the dependent reduced-model accuracy gate.
