# Changes to `src/` made for the SpeedyWeather coupling, and why

Companion to [ecckd_speedyweather.md](ecckd_speedyweather.md). The guiding rule
for that work is to keep changes to the package's core small and to prefer the
extension or SpeedyWeather itself whenever that is possible without giving up
performance. This file lists every change to `src/` on this branch
(`mg/adjust-to-speedy`, Phases 0 and 1 of the plan), what forced it, what the
alternative would have been, and how it is tested. Nothing here changes
numerical results of existing code paths. Later phases (the `EcCKDRadiation`
component and its validation) live on `mg/ecckd-speedy`, whose version of this
file extends the list.

| # | Change | Phase | Lines | Results changed |
|---|---|---|---|---|
| 1 | `AtmosphereProfile`: one array-type parameter per vector | 0 | ~10 | no |
| 2 | `ColumnAtmosphere`: one array-type parameter per array | 1 | ~12 | no |
| 3 | `surface_longwave_emission!` (in place); allocating method delegates | 1 | ~20 | no |
| 4 | `EcCKDTabulatedGasOpticsModel{FT}(model)` element-type conversion | 1 | ~30 | no |
| 5 | Export and API-docs entry for 3, new test file in the runner | 1 | ~3 | no |

## 1. `AtmosphereProfile{NF, VT, VQ, VG}` (was `{NF, V}`)

*File:* `src/column_views.jl`.

*What forced it.* The analytic-band adapter builds an `AtmosphereProfile` from
views into SpeedyWeather arrays. Since SpeedyWeather 0.22 the temperature and
humidity views come from stepped three-dimensional arrays
(`get_prognostic_step(vars.grid.temperature, ...)`) while the geopotential view
comes from a plain two-dimensional field. Their `SubArray` types differ, and
the old struct forced all three vectors to share one type `V`; construction
threw a `MethodError` in `convert`.

*Alternative considered.* Copying the geopotential into a scratch view of the
same shape as the temperature step view. That costs a copy per column per
step for a type-system artefact, and the scratch array would have to mirror
the time-step layout of the prognostic arrays.

*Why acceptable.* Only the type parameters changed; the constructor converts
`surface_pressure`, `rain_rate` and `CO₂` to `NF` as before. No code in the
package dispatched on the second parameter (`radiative_transfer_column.jl` and
`williams_longwave.jl` use `AtmosphereProfile{NF}`).

*Tests.* The extension test `test/test_with_speedyweather.jl` constructs the
profile from real SpeedyWeather views; the solver suite is unchanged.

## 2. `ColumnAtmosphere{FT, PL, PI, TL, TI, G, S, Geo}` (was `{FT, A, G, S, Geo}`)

*File:* `src/runtime_interfaces.jl`.

*What forced it.* The same problem for the staged interface, which the
`EcCKDRadiation` component of the next phase builds per column: layer
temperature is a view into a stepped prognostic array, layer pressure a view
into a `(npoints, nlayers)` work array, interface quantities views into
`(npoints, nlayers + 1)` work arrays. Three different `SubArray` types where
the struct demanded one.

*Alternative considered.* Copying the layer temperature into a work array
every column. Same objection as in 1.

*Why acceptable.* `heating_rates!`, `optical_properties!` and the solvers read
the arrays element-wise and convert to their working precision; none of them
dispatched on the array type. The docstring states that `FT` is the element
type of `temperature_layers`.

*Tests.* `test/test_host_interface.jl`, "ColumnAtmosphere from differently
shaped host views": optics, fluxes and heating rates from views of a 2D/3D host
layout are identical to those from plain vectors.

## 3. `surface_longwave_emission!(out, model, T; emissivity)`

*File:* `src/gas_optics/ecckd_forward.jl`. Exported, documented in
`docs/src/api/ecckd.md`.

*What forced it.* The surface boundary of the tabulated longwave model is
spectral, one value per g-point, and depends on the column's surface
temperature, so a host rebuilds it for every column and (blended over ocean
and land) twice per column. The existing method allocated a `Vector{FT}` per
call. Inside SpeedyWeather's fused column kernel that would be the only
allocation in the hot loop on CPU, and on GPU heap allocation inside a kernel
does not compile at all.

*Alternative considered.* Computing the emission in the extension. That would
require calling the unexported `source_table_bracket` and `longwave_source`,
tying the extension to internals of the source-table interpolation.

*Why acceptable.* A dozen lines; the allocating method now calls the in-place
one, so there is a single implementation. The `f` / `f!` pair is the idiomatic
Julia shape.

*Tests.* `test/test_host_interface.jl`: in-place equals allocating for several
temperatures and emissivities, length check throws, zero allocations on a
Float32 model.

## 4. `EcCKDTabulatedGasOpticsModel{FT}(model)`

*File:* `src/gas_optics/ecckd_forward.jl`.

*What forced it.* SpeedyWeather runs in `Float32` by default, the NetCDF loader
produces `Float64` tables, and `optical_properties!` requires the model and
the optics arrays to share their element type. `Adapt.adapt` moves arrays
between devices but does not change element types.

*Alternative considered.* (a) A `FT` keyword on the NetCDF loader: that lives
in the NCDatasets extension and would still leave already-loaded models
unconvertible. (b) Rebuilding the model field by field in the SpeedyWeather
extension through the keyword constructor: fourteen fields, several optional
with `nothing`/empty-array conventions, i.e. knowledge of the struct layout
that belongs next to the struct.

*Why acceptable.* Thirty lines that only call the existing keyword
constructor, so every validation (log-uniform grids, shapes) is re-run on the
converted arrays. Arrays already of type `FT` are reused. Useful beyond
SpeedyWeather for any single-precision or GPU host.

*Tests.* `test/test_host_interface.jl`: field-wise conversion, reuse for the
same type, absent optional tables stay absent, Float32 optics reproduce
Float64 optics to 1e-4, and reference-sized grids (53-point pressure grid over
five decades, 12-point H2O grid) survive the Float32 re-validation of the
constructor with about a threefold margin on its 1e-5 tolerance.

## 5. Housekeeping

- `src/NumericalRadiation.jl`: exports `surface_longwave_emission!`.
- `docs/src/api/ecckd.md`: `@docs` entry for it (the package builds docs with
  `checkdocs = :exports`).
- `test/runtests.jl`: includes `test_host_interface.jl`.

## Considered and deliberately not changed

- **`@boundscheck` around `check_ecckd_optics_shapes`.** The checks are O(1)
  size comparisons, negligible on CPU, and `@boundscheck` would only remove
  them if `optical_properties!` were inlined into the `@inbounds` caller, which
  it is not. Their `throw` paths need a kernel-safe variant for GPU anyway;
  revisit with the first GPU run.
- **Splitting `optical_properties!` per stream.** Unnecessary once SpeedyWeather
  bundles both streams into one component (upstream change U1): one call fills
  both streams, which is exactly what the component needs.
- **A `(ng, nlayers)`-major work-array layout accessor.** A `PermutedDimsArray`
  of the `(nlayers, ng)` column slice gives the package's `[ig, k]` indexing at
  no cost in the package.
- **Gas-amount and interface-temperature helpers.** Host glue; they belong in
  the extension.
