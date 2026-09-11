# Plan: running ecCKD inside SpeedyWeather.jl

Status legend: `[ ]` open, `[~]` in progress, `[x]` done. Update the checkboxes
as work lands and link PRs next to the items.

Created 2026-09-11. Target versions: NumericalRadiation 0.1.x,
SpeedyWeather 0.22.x plus the upstream changes tracked in
[speedyweather_upstream.md](speedyweather_upstream.md).

Changes that land in SpeedyWeather itself are planned in that companion
file and only referenced here as **U1 … U4**.

Every change made to this package's `src/` for the coupling is listed and
justified in [src_changes_for_speedyweather.md](src_changes_for_speedyweather.md).

## Background

**What exists.** NumericalRadiation already has a complete, host-neutral
ecCKD column path, the staged API in `src/runtime_interfaces.jl`:

1. `read_reference_ecckd_gas_optics` loads an `EcCKDTabulatedGasOpticsModel`
   from the lazy `ecrad_data` artifact (NCDatasets extension).
2. `optical_properties!` (`src/gas_optics/ecckd_forward.jl`) fills
   `LongwaveOptics` and `ShortwaveOptics`, shaped `(ng, nlayers)`, for both
   streams from one pass over the gas tables.
3. `radiative_fluxes!` with `CloudlessLongwave` / `CloudlessShortwave`.
4. `heating_rates!` converts flux convergence to K s⁻¹.

`examples/ecckd_column.jl` is exactly the call sequence the SpeedyWeather
kernel must reproduce per column. Gas inputs are per-layer molar amounts in
mol m⁻², `composite` is the dry-air column, and the H2O mole fraction is
derived from `h2o / composite`.

**What is stale.** `ext/NumericalRadiationSpeedyWeatherExt.jl` only wraps the
analytic-band longwave and targets SpeedyWeather 0.20: it reads
`vars.grid.temperature_prev`, `vars.grid.pressure_prev` and
`model.land_sea_mask.mask`, none of which exist in 0.22. `Project.toml` pins
compat to 0.20.x.

**SpeedyWeather 0.22 contract for a radiation scheme.**

- `parameterization!(ij, vars, scheme, model)` is called once per column,
  fused into one GPU kernel over all parameterizations, in the fixed order
  zenith, albedo, shortwave, longwave.
- Today longwave and shortwave are separate model components; **U1**
  replaces them with one `radiation` component so a scheme can do both
  streams in one call.
- State access: `get_prognostic_step(vars.grid.temperature, model.time_stepping, scheme)`,
  same for humidity; `vars.parameterizations.surface_pressure[ij]`;
  `model.land_sea_mask.land_fraction[ij]`; pressures via `pressure_half` /
  `pressure_thickness(k, pₛ, model.geometry.vertical_coordinates)` (sigma or
  hybrid sigma-pressure).
- Outputs go to `vars.parameterizations.*` fields declared by `variables`.
  Extra per-column work arrays are declared the same way, e.g.
  `Grid4D(n)` allocates `(npoints, nlayers, n)`, `Grid3D(n)` allocates
  `(npoints, n)`; `variables(component, model)` can size them from `model`.
- Available inputs: layer T, specific humidity, pₛ, `cos_zenith`, ocean and
  land albedo, SST, soil temperature, land fraction, scalar CO2 [ppm] in
  `vars.prognostic.greenhouse_gases.co2[]` (if configured), solar constant,
  gravity, heat capacity.
- Not available: ozone (**U2**), CH4, N2O, cloud state, interface
  temperatures, surface emissivity.

**Data.** The lazy `ecrad_data` artifact is not installed locally; the first
run downloads the ecRad source archive (override with `RH_ECRAD_DATA_PATH`).
Which gases the reference files carry must be read from their `gas_names`.
ecCKD climate models normally include h2o, co2, o3, ch4, n2o and CFCs on top
of `composite`. Gases with a reference mole fraction that are omitted from
`names` are implicitly held at that reference value (acceptable default for
CH4 / N2O). Ozone is a linear gas and must be supplied.

## Ordering

```
U1 Radiation bundle (SpeedyWeather, bit-identical)  ──►  release
        │
        ▼
Phase 0 port extension to new SpeedyWeather ──► Phase 1 package prep ──► Phase 2 EcCKDRadiation ──► Phase 4 validation
                                                                                   ▲
U2 ozone, U3 call frequency, U4 array layout (SpeedyWeather, in parallel) ─────────┘
```

Phase 1 does not depend on U1 and can start immediately.

## Phase 0. Port the extension to the post-U1 SpeedyWeather

Done 2026-09-11 against the local U1 branch (SpeedyWeather 0.23.0-DEV, developed
via `Pkg.develop`). Depends on **U1** being *released* only for CI.

- [x] Bump compat in `Project.toml` to `SpeedyWeather = "0.23"`. Pkg accepts the
      developed `0.23.0-DEV` under this compat. Until 0.23 is registered,
      `test/Project.toml` carries a `[sources]` entry pointing at the
      `mg/numericalradiation` branch of SpeedyWeather.jl (monorepo subdir
      `SpeedyWeather`). `[sources]` needs Julia ≥ 1.11, so the package's Julia
      compat and the CI matrix moved from 1.10 to 1.11. Remove the entry once
      SpeedyWeather 0.23 is released.
- [x] Port `SpeedyAnalyticBandLongwave` to 0.22+ accessors
      (`get_prognostic_step` / `get_tendency_step`, `vars.parameterizations.surface_pressure`,
      `land_sea_mask.land_fraction`, `vars.dynamics.geopotential`, stepped SST) and to
      being used as `Radiation(spectral_grid; longwave = SpeedyAnalyticBandLongwave(...))`.
      The kernel is `@propagate_inbounds`.
- [x] Unplanned core change: `AtmosphereProfile` now has one array-type parameter
      per vector (`VT`, `VQ`, `VG`). SpeedyWeather's temperature/humidity views come
      from stepped 3D arrays and the geopotential view from a plain 2D field, so a
      single `V` no longer fits. No other code depended on the two-parameter form.
- [x] `test/test_with_speedyweather.jl` rewritten for the new API (longwave-only
      model, CO₂ forcing, and a 4-step full `run!`); 16 tests pass with
      `--check-bounds=yes`. Core solver (692) and misc (33) tests still pass.
- [x] README "With SpeedyWeather.jl" updated; unused `RingGrids` dep and its
      `=0.1.7` pin dropped from `test/Project.toml`.

## Phase 1. Prepare NumericalRadiation for a fused, allocation-free kernel

Re-scoped 2026-09-11: keep changes to the package small. Only what cannot live
in the extension goes into `src/`; host glue moves to Phase 2. Done on
`mg/adjust-to-speedy`, tests in `test/test_host_interface.jl`.

- [x] ~~Split `optical_properties!` per stream.~~ Not needed: with **U1** a
      single `EcCKDRadiation` component consumes both streams from one call,
      which is exactly what `optical_properties!` produces today.
- [x] In-place `surface_longwave_emission!(out, model, T; emissivity)`
      (exported, allocation-free); the allocating method now calls it.
- [x] Element-type conversion `EcCKDTabulatedGasOpticsModel{FT}(model)`
      (Float64 tables from the loader → Float32 for SpeedyWeather's default
      `NF`; `optical_properties!` requires optics and model to share `FT`).
      Reuses arrays already of type `FT`; absent optional tables stay absent.
- [x] `ColumnAtmosphere` with one array-type parameter per array (as done for
      `AtmosphereProfile` in Phase 0): a host's layer and interface views come
      from arrays of different shape and, for stepped prognostics, different
      rank. No other code depended on the single parameter.
- [ ] ~~Guard `check_ecckd_optics_shapes` behind `@boundscheck`.~~ Deferred:
      the checks are O(1) size comparisons, negligible on CPU, and
      `@boundscheck` would only elide them if `optical_properties!` were
      inlined into the `@inbounds` caller, which it is not. A GPU run needs a
      kernel-safe variant of the `throw` paths anyway; revisit with the first
      GPU test in Phase 4.
- [ ] ~~Layout-agnostic accessor (U4).~~ Deferred to Phase 2: a permuted view
      of the `(nlayers, ng)` column slice keeps the package's `[ig, k]`
      indexing without any package change; only if that costs measurably do
      we revisit.
- [x] Unit tests: in-place emission equals the allocating one and does not
      allocate; converted Float32 model reproduces Float64 optics to 1e-4;
      `ColumnAtmosphere` from views of a 2D/3D host layout gives results
      identical to plain vectors through optics, fluxes and heating rates.

Moved to Phase 2 (extension, host glue): gas amounts from specific humidity,
interface temperatures from layer temperatures.

Re-audit 2026-09-11: no defects found in the committed code. Verified that
reference-sized grids (53-point pressure grid over five decades, 12-point H2O
grid) pass the constructor's 1e-5 log-uniform re-validation after the Float32
round trip with about a threefold margin, and added that as a regression
test. Wrapped the test file in a module like the other consolidated tests, and
corrected the `ColumnAtmosphere` docstring, which overstated that the arrays
"need not share an element type": the kernels do convert on read, but `FT`
remains the working precision. Deferred items unchanged.

## Phase 2. `EcCKDRadiation` SpeedyWeather component in the extension

Implemented 2026-09-11 on `mg/adjust-to-speedy`, clear-sky, CPU.

- [x] `EcCKDRadiation{NF} <: SpeedyWeather.AbstractRadiation` in
      `ext/NumericalRadiationSpeedyWeatherExt/ecckd_radiation.jl` (the extension
      is now a directory: module file, `analytic_band_longwave.jl`,
      `ecckd_radiation.jl`): one component for both
      streams holding the tabulated gas optics (converted to the grid's `NF` at
      construction), `Adapt.@adapt_structure`. Options: `CO₂` default [ppm]
      (used when the model has no `greenhouse_gases.co2`), `ozone` (function of
      pressure or constant; crude Chapman-layer default), `mole_fractions` for
      any further gas of the ecCKD model (required, checked at construction),
      ocean / land emissivity, molar masses. Constructors take a gas-optics
      model, a reference pair name (`"32x32"`, loads via NCDatasets), or nothing
      (default pair). Used as `PrimitiveWetModel(spectral_grid; radiation = EcCKDRadiation(spectral_grid))`.
- [ ] ~~Thin `EcCKDLongwave` / `EcCKDShortwave` wrappers.~~ Not done; low
      priority, they would redo gas optics per stream.
- [x] ~~Load tables in `initialize!`.~~ Changed: tables are loaded and
      converted in the constructor instead. Keeps the struct immutable and
      GPU-adaptable, and `variables` needs the g-point counts before
      `initialize!` runs. NCDatasets stays out of the hot path either way.
- [x] Extension helper `gas_amounts!`: specific humidity, layer Δp, gravity and
      the CO₂ mole fraction → per-layer molar amounts of every gas of the
      model, unrolled over the gas names with a `@generated` function;
      `composite` = dry air, `h2o` from `q`, `co2` from ppm, others from
      `mole_fractions`.
- [x] Extension helper `interface_temperatures!`: linear in pressure between
      layer centres, `T[1]` at the top, the land-fraction-blended surface
      temperature at the bottom (skin temperature, as ecRad/IFS).
- [x] `variables(::EcCKDRadiation, model)`: standard shortwave and longwave
      diagnostics plus a `:ecckd` namespace of work arrays: layer pressure
      (`GridXYZ`), interface pressure and temperature, four interface fluxes,
      two per-g surface-emission vectors (`Grid3D`), gas amounts and the seven
      optical-property arrays (`Grid4D(n = ngas | ng)`), and ten shortwave
      adding-method work arrays. Column views of these build
      `ColumnAtmosphere`, `LongwaveOptics`, `ShortwaveOptics`, `RadiativeFluxes`
      and `CloudlessShortwaveWorkspace` without allocation; the `(nlayers, ng)`
      column slice is wrapped in a `PermutedDimsArray` for `[ig, k]` indexing.
- [x] `parameterization!(ij, vars, ::EcCKDRadiation, model)`, split into
      `ecckd_surface_state`, `ecckd_column_atmosphere!`, `ecckd_column_optics`,
      `ecckd_longwave!`, `ecckd_shortwave!`, `ecckd_heating!`: pressures from
      the vertical coordinates, gas amounts, CO₂ from the greenhouse-gas
      variable, interface temperatures; `optical_properties!` once; longwave
      with spectral surface emission blended over ocean and land (in place, two
      `surface_longwave_emission!` calls); shortwave only for `cos_zenith > 0`
      with TOA down = `solar_constant * cos_zenith` and blended albedo;
      diagnostics including the ocean / land splits; net flux convergence of
      both streams into `dTdt` via SpeedyWeather's `flux_to_tendency`.
- [x] Unplanned package change: the Rayleigh (scattering) path of
      `CloudlessShortwave` allocated twelve vectors per g-point per column, i.e.
      ~400 allocations per column per step, which cannot be avoided from the
      extension without dropping Rayleigh scattering. `radiative_fluxes!` now
      takes an optional `CloudlessShortwaveWorkspace` (ten caller-owned
      vectors; `radiation_workspace(CloudlessShortwave(), optics)` allocates
      one, the keyword constructor accepts host views), and the per-g flux
      accumulation writes directly into the flux arrays instead of through
      scratch vectors. The 5-argument call allocates a workspace only when
      scattering is present. Results are unchanged (checked bitwise in
      `test/test_host_interface.jl`; `test/test_solvers.jl` still passes).
- [x] Clear-sky only in this version. Cloud coupling (package has
      cloud-overlap solvers, SpeedyWeather has no cloud state) is a follow-up.
- [x] Tests in `test/test_with_speedyweather.jl`: construction and work-array
      sizes on the reference 32x32 model, missing-gas error; a realistic column
      (OLR range, surface budget signs, shortwave transmission, moderate heating
      rates, column energy conservation to 1e-3, night = zero shortwave with
      unchanged longwave, 4×CO₂ reduces OLR by a few W/m²); fluxes of one
      column agree with the staged API called directly on the same inputs; a
      4-step full model run.
- [x] README "With SpeedyWeather.jl" updated.

Known limitations of this first version: clear sky; ozone from an analytic
default profile (**U2**); the package's shape checks and broadcasts inside
`optical_properties!` / `radiative_fluxes!` are not yet GPU-safe (Phase 4);
`co2` is a single global value from `greenhouse_gases`.

## Phase 3. Upstream SpeedyWeather changes

Tracked in [speedyweather_upstream.md](speedyweather_upstream.md):

- **U1** one `radiation` component bundling shortwave and longwave,
  bit-identical with the existing schemes first. Prerequisite for Phase 0
  and 2.
- **U2** prescribed ozone.
- **U3** radiation call frequency.
- **U4** work-array layout for g-point arrays (or make NumericalRadiation
  layout-agnostic instead, Phase 1).

## Phase 4. Validation

- [ ] Test: one SpeedyWeather column through `EcCKDRadiation` vs. the same
      inputs through the staged API directly; fluxes agree to tolerance.
- [ ] Aquaplanet run at low resolution with `radiation = EcCKDRadiation(...)`:
      global-mean OLR ≈ 240 W m⁻², closed TOA budget, side-by-side against
      the default `Radiation(OneBandShortwave, OneBandLongwave)`.
- [ ] Per-column benchmark for the 32x32 and 64x96 pairs to size the
      call-frequency decision (**U3**).
- [ ] CI: gate ecCKD tests on the artifact download; reuse
      `RH_ECRAD_DATA_PATH`.

## Open decisions

- Default reference model pair: proposed `climate_32x32`.
- Ozone source until **U2** lands: ship a zonal-mean profile in the extension
  or require the user to pass one.
- Scope of back-compat in **U1** (deprecation shims for the old keywords).
