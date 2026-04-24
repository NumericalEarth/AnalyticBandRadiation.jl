# Single-column comparison: AnalyticBandRadiation (Williams 2026) vs RRTMGP
# clear-sky longwave, on a shared idealized tropical atmosphere.
#
# This script uses Breeze as the host for RRTMGP (via BreezeRRTMGPExt). It
# builds one atmospheric column with Breeze's `AnelasticDynamics` + RRTMGP
# `ClearSkyOptics`, extracts T(z) and qᵥ(z) from the initialized model,
# converts height-centered layers into sigma coordinates, then runs
# `AnalyticBandRadiation.WilliamsLongwave` on the same profile and plots
# both heating-rate profiles together.
#
# Run from the package root with
#     julia --project=examples examples/rrtmgp_comparison.jl
# (the `examples/Project.toml` brings in Breeze, RRTMGP, NCDatasets, etc.)

import AnalyticBandRadiation as ABR
using Breeze
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znodes
using NCDatasets
using RRTMGP
using ClimaComms
using CairoMakie
using Dates
using Statistics

# ---------------------------------------------------------------------------
# 1. Build a single-column Breeze model with RRTMGP clear-sky longwave.
# ---------------------------------------------------------------------------

Nz = 64
λ, φ = -76.13, 39.48           # Maryland; irrelevant here but needed for zenith
grid = RectilinearGrid(size = Nz, x = λ, y = φ, z = (0, 20kilometers),
                       topology = (Flat, Flat, Bounded))

surface_temperature = 300.0
constants = ThermodynamicConstants()
reference_state = ReferenceState(grid, constants;
                                 surface_pressure = 101_325,
                                 potential_temperature = surface_temperature)
dynamics = AnelasticDynamics(reference_state)
microphysics = SaturationAdjustment(equilibrium = WarmPhaseEquilibrium())

# Clear-sky full-spectrum RRTMGP reference.
radiation_rrtmgp = RadiativeTransferModel(grid, ClearSkyOptics(), constants;
                                           surface_temperature,
                                           surface_emissivity = 0.98,
                                           surface_albedo = 0.1,
                                           solar_constant = 1361.0)

clock = Clock(time = DateTime(1950, 11, 1, 12, 0, 0))
model = AtmosphereModel(grid; clock, dynamics, microphysics, radiation = radiation_rrtmgp)

# Tropical moisture profile.
θ₀   = reference_state.potential_temperature
q₀   = 0.020
Hᵗ   = 3000.0
qᵗᵢ(z) = q₀ * exp(-z / Hᵗ)
set!(model; θ = θ₀, qᵗ = qᵗᵢ)

# ---------------------------------------------------------------------------
# 2. Extract the column state into plain vectors.
# ---------------------------------------------------------------------------

# RRTMGP flux divergence → heating rate via -dF/dz / (ρ cₚ). Breeze stores the
# pre-computed flux_divergence on cell centers in W/m³.
ρ_ref = Array(interior(reference_state.density, 1, 1, :))
cₚ    = constants.dry_air.heat_capacity
flux_div_rrtmgp = Array(interior(radiation_rrtmgp.flux_divergence, 1, 1, :))
heating_rate_rrtmgp = flux_div_rrtmgp ./ (ρ_ref .* cₚ) .* 86400   # K/day

# Temperature and humidity at cell centers (height coordinate).
T      = Array(interior(model.temperature, 1, 1, :))
qᵛ_field = specific_humidity(model)
qᵛ     = Array(interior(qᵛ_field, 1, 1, :))
p_ref  = Array(interior(reference_state.pressure, 1, 1, :))
z_full = znodes(grid, Center())

# Ordering convention: AnalyticBandRadiation uses top-down (k=1 is TOA). Breeze
# uses bottom-up (k=1 is surface). Reverse now.
T_top_down      = reverse(T)
qᵛ_top_down     = reverse(qᵛ)
p_full_top_down = reverse(p_ref)

surface_pressure = p_ref[1]   # at z = 0 (bottom of domain in Breeze ordering)

# ---------------------------------------------------------------------------
# 3. Build sigma coordinates from the reference-state pressure profile.
# ---------------------------------------------------------------------------

# Half-level pressures: linear interpolation of p(z) to the cell faces.
p_half = similar(p_ref, Nz + 1)
p_half[1]   = surface_pressure                 # z = 0
p_half[end] = 2 * p_ref[end] - p_ref[end - 1]  # extrapolate to model top
for k in 2:Nz
    p_half[k] = (p_ref[k - 1] + p_ref[k]) / 2
end

σ_half = reverse(p_half) ./ surface_pressure   # top-down, monotonically 0 → 1
geometry = ABR.ColumnGeometry(σ_half)

profile = ABR.ColumnProfile(
    temperature      = T_top_down,
    humidity         = qᵛ_top_down,
    geopotential     = zeros(Nz),
    surface_pressure = surface_pressure,
)

surface = ABR.ColumnSurface{Float64}(sea_surface_temperature  = surface_temperature,
                                      land_surface_temperature = NaN,
                                      land_fraction            = 0.0)
pc = ABR.PhysicalConstants{Float64}(heat_capacity = cₚ)

lw_williams = ABR.WilliamsLongwave(Float64; do_co2 = true, co2_ppmv = 420.0)
dT_williams = zeros(Nz)
diag = ABR.LongwaveDiagnostics{Float64}()
ABR.solve_longwave!(dT_williams, diag, lw_williams, profile, geometry, surface, pc)

heating_rate_williams = reverse(dT_williams) .* 86400   # flip back to bottom-up

# ---------------------------------------------------------------------------
# 4. Plot the comparison.
# ---------------------------------------------------------------------------

fig = Figure(size = (900, 500))
ax1 = Axis(fig[1, 1]; xlabel = "Temperature (K)", ylabel = "Altitude (km)",
           title = "Shared atmospheric profile")
lines!(ax1, T, z_full ./ 1000.0; color = :gray30, linewidth = 2)

ax2 = Axis(fig[1, 2]; xlabel = "Clear-sky LW heating rate (K day⁻¹)",
           ylabel = "Altitude (km)", title = "Williams (2026) vs RRTMGP")
lines!(ax2, heating_rate_rrtmgp,  z_full ./ 1000.0;
       color = :dodgerblue, linewidth = 2.5, label = "RRTMGP clear-sky")
lines!(ax2, heating_rate_williams, z_full ./ 1000.0;
       color = :crimson,    linewidth = 2.5, label = "Williams (2026) SSM")
vlines!(ax2, 0; color = :gray70, linestyle = :dash)
axislegend(ax2; position = :rt)

save(joinpath(@__DIR__, "rrtmgp_comparison.png"), fig)

# ---------------------------------------------------------------------------
# 5. Scalar diagnostics.
# ---------------------------------------------------------------------------

# TOA outgoing LW from RRTMGP = upwelling flux at the model top.
olr_rrtmgp  = Array(interior(radiation_rrtmgp.upwelling_longwave_flux, 1, 1, :))[end]
olr_williams = diag.outgoing_longwave

# Vertically integrated heating rate (K/day, mass-weighted mean).
Δz_c = diff(znodes(grid, Face()))
mass = ρ_ref .* Δz_c
mean_hr_rrtmgp  = sum(heating_rate_rrtmgp  .* mass) / sum(mass)
mean_hr_williams = sum(heating_rate_williams .* mass) / sum(mass)

println("OLR (RRTMGP)        = ", round(olr_rrtmgp;  digits = 2), " W/m²")
println("OLR (Williams SSM) = ", round(olr_williams; digits = 2), " W/m²")
println("Mean column heating rate:")
println("  RRTMGP  = ", round(mean_hr_rrtmgp;  digits = 3), " K/day")
println("  Williams = ", round(mean_hr_williams; digits = 3), " K/day")
