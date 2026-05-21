# Architecture

The package exposes the existing analytic-band column solvers through three
levels of API.

## High-level column update

[`radiative_heating!`](@ref) runs the reusable [`RadiativeTransferColumn`](@ref)
workspace by calling longwave and shortwave component solvers. This is the
convenience path for examples and single-column workflows.

## Component access

Host models can call [`solve_longwave!`](@ref), [`solve_shortwave!`](@ref), and
[`heating_rates!`](@ref) independently. This keeps solver, flux, and tendency
ownership explicit for models that need their own vertical integrals or
tendency insertion.

The staged interface also defines [`optical_properties!`](@ref),
[`cloud_optical_properties!`](@ref), [`aerosol_optical_properties!`](@ref), and
[`radiative_fluxes!`](@ref) for future gas-optics and solver implementations.
The first concrete staged solver is [`CloudlessLongwave`](@ref), which consumes
[`LongwaveOpticalProperties`](@ref) and writes caller-owned
[`RadiativeFluxes`](@ref). [`CloudlessShortwave`](@ref) does the same for
absorptive shortwave optical depths stored in [`ShortwaveOpticalProperties`](@ref).

[`heating_rates!`](@ref) can convert [`RadiativeFluxes`](@ref) into layer
heating rates for a [`ColumnAtmosphere`](@ref) using explicit `gravity` and
`heat_capacity` keywords. The convention is top-down pressure interfaces,
net-downward flux, and positive heating for atmospheric warming.

## Workspace access

[`radiation_workspace`](@ref) returns reusable storage for repeated runtime
calls. For the current analytic-band path, [`RadiativeTransferColumn`](@ref) is
already the workspace: it owns the temperature-tendency vector, shortwave
transmissivity scratch, and diagnostics.

Future Breeze and SpeedyWeather integrations should prefer caller-owned arrays,
views, or explicit workspaces over per-call temporary arrays.
