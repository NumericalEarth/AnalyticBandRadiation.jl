# Notation

`AnalyticBandRadiation` follows the [NumericalEarth.jl notation
guide](https://github.com/NumericalEarth/NumericalEarth.jl/blob/main/docs/src/appendix/notation.md)
for symbolic names in math, in docstring equations, and in plot labels. Julia
struct fields stay descriptive snake_case (matching Breeze's public API
convention), but the mapping to the symbolic notation is always unambiguous.

## Radiative fluxes

| Math | Where it appears in code | NumericalEarth symbol |
|:-----|:-------------------------|:----------------------|
| `ℐꜜˡʷ` at surface | [`LongwaveDiagnostics`](@ref) `.surface_longwave_down` | `\scrI\^downarrow\^l\^w` |
| `ℐꜛˡʷ` at surface | [`LongwaveDiagnostics`](@ref) `.surface_longwave_up`, `.ocean_surface_longwave_up`, `.land_surface_longwave_up` | `\scrI\^uparrow\^l\^w` |
| `ℐꜛˡʷ` at TOA (OLR) | [`LongwaveDiagnostics`](@ref) `.outgoing_longwave` | `\scrI\^uparrow\^l\^w` at top of atmosphere |
| `ℐꜜˢʷ` at surface | [`ShortwaveDiagnostics`](@ref) `.surface_shortwave_down` | `\scrI\^downarrow\^s\^w` |
| `ℐꜛˢʷ` at surface | [`ShortwaveDiagnostics`](@ref) `.surface_shortwave_up` | `\scrI\^uparrow\^s\^w` |
| `ℐꜛˢʷ` at TOA | [`ShortwaveDiagnostics`](@ref) `.outgoing_shortwave` | reflected-to-space shortwave |

## State and surface variables

| Math | Code | Description |
|:-----|:-----|:------------|
| `T` | `ColumnProfile.temperature` | Layer air temperature (K) |
| `q` | `ColumnProfile.humidity` | Specific humidity (kg kg⁻¹) |
| `p` | `ColumnProfile.surface_pressure`; `ColumnGeometry.σ_*·pₛ` in-code | Pressure (Pa) |
| `α` | `ColumnSurface.ocean_albedo`, `land_albedo` | Surface albedo |
| `ϵ` | `ColumnSurface.ocean_emissivity`, `land_emissivity` | Surface emissivity |
| `σ` (Stefan–Boltzmann) | `PhysicalConstants.stefan_boltzmann` | W m⁻² K⁻⁴ |
| `g` | `PhysicalConstants.gravity` | m s⁻² |
| `cₚ` | `PhysicalConstants.heat_capacity` | Isobaric specific heat (J kg⁻¹ K⁻¹) |

## Sigma-coordinate vertical grid

The vertical grid is stored in [`ColumnGeometry`](@ref) using the
pressure-normalized `σ = p / pₛ` convention inherited from SpeedyWeather.

| Math | Code | Description |
|:-----|:-----|:------------|
| `σₖ` at midpoints | `σ_full` | Length `nlayers` |
| `σₖ₊½` at interfaces | `σ_half` | Length `nlayers + 1`, monotonic from 0 (TOA) to 1 (surface) |
| `Δσₖ` | `σ_thick` | Layer thickness, `= diff(σ_half)` |

Note: the package uses `σ` for the vertical coordinate *and* `σ` in context
for the Stefan–Boltzmann constant (e.g. `PhysicalConstants.stefan_boltzmann`).
Local variables in the solvers disambiguate with `σ_SB`.

## Longwave spectroscopy

| Math | Code | Description |
|:-----|:-----|:------------|
| `ν̃` | `ν̃` | Wavenumber (cm⁻¹) |
| `B(T, ν̃)` | [`planck_wavenumber`](@ref) | Spectral Planck radiance |
| `κ_line^ref(ν̃)` | [`h2o_line_kappa_ref`](@ref) | Reference H₂O line absorption |
| `κ_cnt^ref(ν̃)` | [`h2o_cont_kappa_ref`](@ref) | Reference H₂O continuum absorption |
| `κ_CO₂^ref(ν̃)` | [`co2_kappa_ref`](@ref) | Reference CO₂ absorption |
| `τ(p)` | integrated internally | Optical depth |
| `D` | `WilliamsLongwave.diffusivity` | Two-stream diffusivity factor (≈ 1.5) |

Loop-internal scalar accumulators in [`solve_longwave!`](@ref) use compact
names `U` and `D` for `ℐꜛ` and `ℐꜜ` in the spectral sweep. Comments in the
solver tie them to the math.
