# Shortwave: SPEEDY one-band scheme

The [`OneBandShortwave`](@ref) solver composes three sub-schemes — a
diagnostic cloud model, a layer transmissivity, and a radiative-transfer
solver — into a single call. It reproduces Fortran SPEEDY (Kucharski,
Molteni & Bracco, 2006, Appendix B) with the same parameter defaults.

## Transmissivity sensitivity to zenith angle

```@example sw_zenith
using AnalyticBandRadiation
using CairoMakie

nlayers = 32
σ_half  = collect(range(0.0, 1.0, length = nlayers + 1))
geom    = ColumnGrid(σ_half)

profile = AtmosphereProfile(
    temperature      = collect(range(220.0, 295.0, length = nlayers)),
    humidity         = fill(0.005, nlayers),
    geopotential     = zeros(nlayers),
    surface_pressure = 100_000.0,
)
constants = PhysicalConstants{Float64}()
thermo    = ThermodynamicConstants{Float64}()
scheme    = OneBandShortwave(Float64)

cos_zeniths = [0.2, 0.4, 0.6, 0.8, 1.0]
fig = Figure(size = (780, 420))
ax  = Axis(fig[1, 1];
           xlabel = "SW heating rate [K day⁻¹]",
           ylabel = "σ (nondimensional pressure)",
           yreversed = true,
           title  = "Shortwave heating rate vs zenith angle")

for μ in cos_zeniths
    surface = SurfaceState{Float64}(sea_surface_temperature = 295.0,
                                     land_surface_temperature = NaN,
                                     land_fraction = 0.0,
                                     ocean_albedo = 0.07,
                                     land_albedo  = 0.07,
                                     cos_zenith   = μ)
    dT   = zeros(nlayers)
    dg   = ShortwaveDiagnostics{Float64}(nlayers)
    tbuf = similar(profile.temperature)
    solve_shortwave!(dT, dg, scheme, profile, geom, surface, constants, thermo;
                     transmissivity_scratch = tbuf)
    lines!(ax, dT .* 86400, geom.σ_full; label = "cos θ = $μ", linewidth = 2)
end
axislegend(ax; position = :rb)
save("sw_zenith.png", fig); nothing # hide
```

![](sw_zenith.png)

Near the TOA the heating is dominated by ozone; near the surface by water
vapour. At grazing angles (cos θ → 0) the optical path is long and heating
is concentrated aloft, matching the SPEEDY zenith-correction factor
`(1 + a_zen (1 − cos θ)^n_zen)` in [`BackgroundShortwaveTransmissivity`](@ref).

## Cloud-albedo sensitivity

Repeating the same column but sweeping the cloud cover through the
diagnostic scheme's precipitation term shows the surface-insolation response.

```@example sw_clouds
using AnalyticBandRadiation
using CairoMakie

nlayers = 16
σ_half  = collect(range(0.0, 1.0, length = nlayers + 1))
geom    = ColumnGrid(σ_half)
base_profile = AtmosphereProfile(
    temperature      = collect(range(220.0, 295.0, length = nlayers)),
    humidity         = fill(0.008, nlayers),
    geopotential     = zeros(nlayers),
    surface_pressure = 100_000.0,
)
surface = SurfaceState{Float64}(sea_surface_temperature = 295.0,
                                 land_surface_temperature = NaN,
                                 land_fraction = 0.0,
                                 ocean_albedo = 0.07,
                                 land_albedo  = 0.07,
                                 cos_zenith   = 0.6)
constants = PhysicalConstants{Float64}()
thermo    = ThermodynamicConstants{Float64}()

rain_rates = [0.0, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4]   # m/s
ssr   = Float64[]
olw_r = Float64[]
covers = Float64[]
for r in rain_rates
    profile = AtmosphereProfile(temperature = base_profile.temperature,
                            humidity    = base_profile.humidity,
                            geopotential = base_profile.geopotential,
                            surface_pressure = base_profile.surface_pressure,
                            rain_rate = r)
    dT   = zeros(nlayers)
    dg   = ShortwaveDiagnostics{Float64}(nlayers)
    tbuf = similar(profile.temperature)
    solve_shortwave!(dT, dg, OneBandShortwave(Float64), profile, geom,
                     surface, constants, thermo; transmissivity_scratch = tbuf)
    push!(ssr, dg.surface_shortwave_down)
    push!(olw_r, dg.outgoing_shortwave)
    push!(covers, dg.cloud_cover)
end

fig = Figure(size = (780, 360))
ax  = Axis(fig[1, 1];
           xlabel = "Diagnosed cloud cover",
           ylabel = "Flux [W m⁻²]",
           title  = "Cloud-cover response of SW fluxes")
scatter!(ax, covers, ssr;   label = "ℐꜜˢʷ at surface",    markersize = 10)
scatter!(ax, covers, olw_r; label = "ℐꜛˢʷ at TOA",         markersize = 10)
axislegend(ax; position = :rb)
save("sw_clouds.png", fig); nothing # hide
```

![](sw_clouds.png)

As cloud cover grows, more solar flux is reflected back to space (TOA up
rises) and less reaches the surface.
