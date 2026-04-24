# AnalyticBandRadiation.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://NumericalEarth.github.io/AnalyticBandRadiation.jl/dev/)
[![CI](https://github.com/NumericalEarth/AnalyticBandRadiation.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NumericalEarth/AnalyticBandRadiation.jl/actions/workflows/CI.yml)
[![Documenter](https://github.com/NumericalEarth/AnalyticBandRadiation.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/NumericalEarth/AnalyticBandRadiation.jl/actions/workflows/Documenter.yml)

Analytic-band atmospheric radiation for intermediate-complexity models.

This package bundles two per-column radiation schemes:

- **Longwave** — Williams (2026) *Simple Spectral Model*: a 41-wavenumber
  clear-sky two-stream Schwarzschild solver with analytic H₂O line,
  H₂O continuum and CO₂ absorption coefficients.
  Published in *J. Adv. Model. Earth Syst.*, doi:[10.1029/2025MS005405](https://doi.org/10.1029/2025MS005405).
- **Shortwave** — a one-band scheme after SPEEDY (Kucharski, Molteni & Bracco,
  *Quart. J. Roy. Meteor. Soc.*, 2006), with transparent, constant, or
  Kucharski-style background-transmissivity options and a diagnostic
  cloud/stratocumulus model.

The column solvers are pure scalar operations — no allocations, no host-side
loops, GPU-safe — so the same code can be driven by
[SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) (which
dispatches per column via `parameterization!(ij, vars, scheme, model)`) or by
[Breeze](https://github.com/NumericalEarth/Breeze.jl) (which launches a
KernelAbstractions kernel per column inside `_update_radiation!`).

Package extensions wire each host automatically when both packages are loaded.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/NumericalEarth/AnalyticBandRadiation.jl")
```

## Standalone usage (single column)

```julia
using AnalyticBandRadiation

nlayers = 8
σ_half  = collect(range(0.0, 1.0, length = nlayers + 1))
geom    = ColumnGeometry(σ_half)

# Lapse-rate profile: top of atmosphere (k=1) cold, surface (k=nlayers) warm.
profile = ColumnProfile(
    temperature      = collect(range(220.0, 295.0, length = nlayers)),
    humidity         = fill(0.005, nlayers),
    geopotential     = zeros(nlayers),
    surface_pressure = 100_000.0,
)

surface   = ColumnSurface{Float64}(
    sea_surface_temperature  = 295.0,
    land_surface_temperature = 285.0,
    land_fraction            = 0.3,
    ocean_albedo             = 0.07,
    land_albedo              = 0.25,
    cos_zenith               = 0.5,
)
constants = PhysicalConstants{Float64}()

# Longwave
lw   = WilliamsLongwave(Float64; do_co2 = true, co2_ppmv = 280.0)
dTdt = zeros(nlayers)
diag = LongwaveDiagnostics{Float64}()
solve_longwave!(dTdt, diag, lw, profile, geom, surface, constants)

@show diag.outgoing_longwave      # W m⁻²
@show diag.surface_longwave_down  # W m⁻²

# Shortwave
sw    = OneBandShortwave(Float64)
dTdt  = zeros(nlayers)
sdiag = ShortwaveDiagnostics{Float64}(nlayers)
tbuf  = similar(profile.temperature)
thermo = ThermodynamicConstants{Float64}()
solve_shortwave!(dTdt, sdiag, sw, profile, geom, surface, constants, thermo;
                 transmissivity_scratch = tbuf)
```

## With SpeedyWeather.jl

The `AnalyticBandRadiationSpeedyWeatherExt` extension defines a
`SpeedyWilliamsLongwave` scheme that subtypes `SpeedyWeather.AbstractLongwave`
and can be passed directly to `PrimitiveWetModel`:

```julia
using SpeedyWeather, AnalyticBandRadiation
const SpeedyExt = Base.get_extension(AnalyticBandRadiation,
                                     :AnalyticBandRadiationSpeedyWeatherExt)

spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
longwave      = SpeedyExt.SpeedyWilliamsLongwave(spectral_grid; do_co2 = true)
model         = PrimitiveWetModel(spectral_grid; longwave_radiation = longwave)
```

The extension also exports `SimpleSpectralLongwave` as an alias of
`SpeedyWilliamsLongwave` for drop-in compatibility with SpeedyWeather PR #1057.

## With Breeze

The Breeze extension is a scaffold: it defines `SpectralOptics`, the
`RadiativeTransferModel(grid, ::SpectralOptics, constants; …)` constructor,
and the kernel skeleton. Full integration (ZFaceField flux allocation,
reference-pressure plumbing, surface albedo handling) is tracked as a
follow-up and will mirror the
[BreezeRRTMGPExt](https://github.com/NumericalEarth/Breeze.jl/tree/main/ext/BreezeRRTMGPExt)
gray-radiation constructor.

## Schemes at a glance

| Scheme | Purpose | References |
|---|---|---|
| `WilliamsLongwave` | 41-band clear-sky LW | Williams (2026); Armstrong (1968); Mlawer et al. (1997) |
| `TransparentShortwave` | Zero-atmosphere SW; surface-albedo only | — |
| `OneBandShortwave` | SPEEDY moist SW (diagnostic clouds + background transmissivity) | Kucharski, Molteni & Bracco (2006) |
| `OneBandGreyShortwave` | SPEEDY dry SW (no clouds, constant transmissivity) | Kucharski, Molteni & Bracco (2006) |
| `DiagnosticClouds` | Cloud cover from RH + precipitation; stratocumulus from DSE stability | SPEEDY §B4 |
| `BackgroundShortwaveTransmissivity` | Dry-air + aerosol + WV + cloud absorptivities, pressure-weighted | SPEEDY §B4 |
| `ConstantShortwaveTransmissivity` | Single-value column transmissivity | — |

## Tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

## License

MIT. See [LICENSE](./LICENSE).
