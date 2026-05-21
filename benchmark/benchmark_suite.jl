using AnalyticBandRadiation
using Dates
using Printf
using Statistics

function analytic_column(nlayers::Integer = 32)
    grid = ColumnGrid(collect(range(0.0, 1.0, length = nlayers + 1)))
    profile = AtmosphereProfile(
        temperature = collect(range(215.0, 295.0, length = nlayers)),
        humidity = [0.018 * exp(-z / 2.0) for z in range(0.0, 8.0, length = nlayers)],
        geopotential = zeros(nlayers),
        surface_pressure = 100_000.0,
        CO₂ = 420.0,
    )
    surface = SurfaceState(
        sea_surface_temperature = 295.0,
        land_surface_temperature = 288.0,
        land_fraction = 0.25,
        ocean_albedo = 0.07,
        land_albedo = 0.22,
        cos_zenith = 0.5,
    )
    return RadiativeTransferColumn(; grid, profile, surface)
end

function percentile(sorted_values, p)
    i = clamp(ceil(Int, p * length(sorted_values)), 1, length(sorted_values))
    return sorted_values[i]
end

function run_case(; samples = 25, nlayers = 32)
    rtm = analytic_column(nlayers)
    radiative_heating!(rtm)

    times_ms = Float64[]
    for _ in 1:samples
        elapsed = @elapsed radiative_heating!(rtm)
        push!(times_ms, 1e3 * elapsed)
    end

    sorted = sort(times_ms)
    allocations = @allocated radiative_heating!(rtm)

    return (
        case = "analytic_column",
        date = string(Dates.now()),
        backend = "CPU",
        precision = string(eltype(rtm.temperature_tendency)),
        columns = 1,
        layers = nlayers,
        runtime_ms_minimum = minimum(times_ms),
        runtime_ms_median = median(times_ms),
        runtime_ms_p90 = percentile(sorted, 0.90),
        allocations = allocations,
    )
end

function staged_ecckd_column(nlayers::Integer = 32)
    pressure_interfaces = collect(range(1_000.0, 100_000.0, length = nlayers + 1))
    pressure_layers = @views (pressure_interfaces[1:end-1] .+ pressure_interfaces[2:end]) ./ 2
    temperature_interfaces = collect(range(210.0, 298.0, length = nlayers + 1))
    temperature_layers = @views (temperature_interfaces[1:end-1] .+ temperature_interfaces[2:end]) ./ 2
    height = collect(range(0.0, 1.0, length = nlayers))

    atmosphere = ColumnAtmosphere(
        pressure_layers = pressure_layers,
        pressure_interfaces = pressure_interfaces,
        temperature_layers = temperature_layers,
        temperature_interfaces = temperature_interfaces,
        gases = (
            h2o = [0.018 * exp(-8z) for z in height],
            co2 = 420.0e-6,
            o3 = [8e-6 * exp(-((z - 0.35) / 0.18)^2) for z in height],
        ),
        surface = (; temperature = 298.0, emissivity = 1.0, albedo = 0.08),
        geometry = (; cos_zenith = 0.5),
    )

    gas_model = EcCKDGasOpticsModel(
        gas_names = (:h2o, :co2, :o3),
        longwave_absorption = [7.0 0.18 0.70;
                               3.0 0.10 1.20;
                               1.0 0.05 0.30;
                               0.4 0.02 0.10],
        shortwave_absorption = [0.16 0.005 0.20;
                                0.07 0.002 0.45;
                                0.02 0.001 0.30],
        longwave_source_scale = [0.90, 1.00, 1.05, 1.10],
        longwave_weights = [0.20, 0.30, 0.30, 0.20],
        shortwave_weights = [0.35, 0.40, 0.25],
    )

    longwave = LongwaveOpticalProperties(zeros(4, nlayers), zeros(4, nlayers);
                                         weights = zeros(4))
    shortwave = ShortwaveOpticalProperties(zeros(3, nlayers); weights = zeros(3))
    fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    heating = zeros(nlayers)
    lw_boundary = LongwaveBoundaryConditions(
        surface_longwave_up = 5.670374419e-8 * atmosphere.surface.temperature^4,
    )
    sw_boundary = ShortwaveBoundaryConditions(
        toa_shortwave_down = 1361.0 * atmosphere.geometry.cos_zenith,
        surface_albedo = atmosphere.surface.albedo,
    )

    return (; atmosphere, gas_model, longwave, shortwave, fluxes, heating,
            lw_boundary, sw_boundary)
end

function median_runtime_ms!(f, samples)
    f()
    times = Float64[]
    for _ in 1:samples
        push!(times, 1e3 * @elapsed f())
    end
    return median(times), percentile(sort(times), 0.90)
end

function run_staged_ecckd_case(; samples = 25, nlayers = 32)
    state = staged_ecckd_column(nlayers)

    gas!() = optical_properties!(state.longwave, state.shortwave,
                                 state.gas_model, state.atmosphere)
    lw!() = radiative_fluxes!(state.fluxes, CloudlessLongwave(), state.longwave,
                              state.atmosphere, state.lw_boundary)
    sw!() = radiative_fluxes!(state.fluxes, CloudlessShortwave(), state.shortwave,
                              state.atmosphere, state.sw_boundary)
    heat!() = heating_rates!(state.heating, state.fluxes, state.atmosphere;
                             gravity = 9.80665, heat_capacity = 1004.0)
    pipeline!() = begin
        gas!()
        lw!()
        sw!()
        heat!()
    end

    pipeline!()
    gas_ms, gas_p90 = median_runtime_ms!(gas!, samples)
    lw_ms, lw_p90 = median_runtime_ms!(lw!, samples)
    sw_ms, sw_p90 = median_runtime_ms!(sw!, samples)
    heating_ms, heating_p90 = median_runtime_ms!(heat!, samples)
    total_ms, total_p90 = median_runtime_ms!(pipeline!, samples)

    return (
        case = "staged_ecckd_cloudless_column",
        date = string(Dates.now()),
        backend = "CPU",
        precision = string(eltype(state.heating)),
        columns = 1,
        layers = nlayers,
        lw_gpoints = size(state.longwave.optical_depth, 1),
        sw_gpoints = size(state.shortwave.optical_depth, 1),
        gas_optics_runtime_ms_median = gas_ms,
        gas_optics_runtime_ms_p90 = gas_p90,
        longwave_solver_runtime_ms_median = lw_ms,
        longwave_solver_runtime_ms_p90 = lw_p90,
        shortwave_solver_runtime_ms_median = sw_ms,
        shortwave_solver_runtime_ms_p90 = sw_p90,
        heating_runtime_ms_median = heating_ms,
        heating_runtime_ms_p90 = heating_p90,
        total_runtime_ms_median = total_ms,
        total_runtime_ms_p90 = total_p90,
        gas_optics_allocations = (@allocated gas!()),
        longwave_solver_allocations = (@allocated lw!()),
        shortwave_solver_allocations = (@allocated sw!()),
        heating_allocations = (@allocated heat!()),
        total_allocations = (@allocated pipeline!()),
    )
end

function rcemip_style_gas_model()
    pressure_grid = collect(range(1_000.0, 100_000.0, length = 4))
    temperature_grid = collect(range(200.0, 310.0, length = 4))
    longwave_table = zeros(4, 3, length(pressure_grid), length(temperature_grid))
    shortwave_table = zeros(3, 3, length(pressure_grid), length(temperature_grid))

    for ig in axes(longwave_table, 1), j in axes(longwave_table, 2),
        ip in axes(longwave_table, 3), it in axes(longwave_table, 4)
        longwave_table[ig, j, ip, it] =
            0.08 / ig + 0.01j + 2.0e-7pressure_grid[ip] +
            3.0e-5temperature_grid[it]
    end
    for ig in axes(shortwave_table, 1), j in axes(shortwave_table, 2),
        ip in axes(shortwave_table, 3), it in axes(shortwave_table, 4)
        shortwave_table[ig, j, ip, it] =
            0.03 / ig + 0.004j + 1.0e-7pressure_grid[ip] +
            2.0e-5temperature_grid[it]
    end

    return EcCKDTabulatedGasOpticsModel(
        gas_names = (:h2o, :co2, :o3),
        pressure_grid = pressure_grid,
        temperature_grid = temperature_grid,
        longwave_absorption = longwave_table,
        shortwave_absorption = shortwave_table,
        longwave_source_scale = [0.9, 1.0, 1.05, 1.1],
        longwave_weights = [0.20, 0.30, 0.30, 0.20],
        shortwave_weights = [0.35, 0.40, 0.25],
    )
end

function rcemip_style_state(ix, iy, nx, ny, nlayers, gas_model)
    η_interface = collect(range(0.0, 1.0, length = nlayers + 1))
    η_layer = @views (η_interface[1:end-1] .+ η_interface[2:end]) ./ 2
    x = (ix - 1) / max(nx - 1, 1)
    y = (iy - 1) / max(ny - 1, 1)
    surface_temperature = 298.0 + 2.0sin(2π * x) * cos(2π * y)
    surface_pressure = 100_000.0 - 600.0cos(2π * x) * sin(2π * y)

    pressure_interfaces = [1_000.0 + (surface_pressure - 1_000.0) * η^1.15
                           for η in η_interface]
    pressure_layers = @views (pressure_interfaces[1:end-1] .+ pressure_interfaces[2:end]) ./ 2
    temperature_layers = [205.0 + (surface_temperature - 205.0) * η^0.55 -
                          4.0sin(π * η) * cos(2π * x) for η in η_layer]
    temperature_interfaces = [205.0 + (surface_temperature - 205.0) * η^0.55
                              for η in η_interface]
    h2o = [0.020 * exp(-7η) * (1 + 0.15sin(2π * x) * cos(2π * y))
           for η in η_layer]
    o3 = [8.0e-6 * exp(-((η - 0.35) / 0.18)^2) for η in η_layer]

    atmosphere = ColumnAtmosphere(
        pressure_layers = pressure_layers,
        pressure_interfaces = pressure_interfaces,
        temperature_layers = temperature_layers,
        temperature_interfaces = temperature_interfaces,
        gases = (; h2o, co2 = 420.0e-6, o3),
        surface = (; temperature = surface_temperature, emissivity = 1.0, albedo = 0.07),
        geometry = (; cos_zenith = 0.5),
    )
    longwave = LongwaveOpticalProperties(zeros(4, nlayers), zeros(4, nlayers);
                                         weights = zeros(4))
    shortwave = ShortwaveOpticalProperties(zeros(3, nlayers); weights = zeros(3))
    fluxes = RadiativeFluxes(
        longwave_up = zeros(nlayers + 1),
        longwave_down = zeros(nlayers + 1),
        shortwave_up = zeros(nlayers + 1),
        shortwave_down = zeros(nlayers + 1),
    )
    heating = zeros(nlayers)
    lw_boundary = LongwaveBoundaryConditions(
        surface_longwave_up = 5.670374419e-8 * surface_temperature^4,
    )
    sw_boundary = ShortwaveBoundaryConditions(
        toa_shortwave_down = 1361.0 * atmosphere.geometry.cos_zenith,
        surface_albedo = atmosphere.surface.albedo,
    )
    return (; atmosphere, gas_model, longwave, shortwave, fluxes, heating,
            lw_boundary, sw_boundary)
end

function rcemip_style_workload(; nx = 16, ny = 16, nlayers = 64)
    gas_model = rcemip_style_gas_model()
    states = [rcemip_style_state(ix, iy, nx, ny, nlayers, gas_model)
              for iy in 1:ny for ix in 1:nx]
    return (; states, nx, ny, nlayers)
end

function run_rcemip_style_case(; samples = 7, nx = 16, ny = 16, nlayers = 64)
    workload = rcemip_style_workload(; nx, ny, nlayers)
    states = workload.states

    gas!() = foreach(states) do state
        optical_properties!(state.longwave, state.shortwave, state.gas_model, state.atmosphere)
    end
    lw!() = foreach(states) do state
        radiative_fluxes!(state.fluxes, CloudlessLongwave(), state.longwave,
                          state.atmosphere, state.lw_boundary)
    end
    sw!() = foreach(states) do state
        radiative_fluxes!(state.fluxes, CloudlessShortwave(), state.shortwave,
                          state.atmosphere, state.sw_boundary)
    end
    heat!() = foreach(states) do state
        heating_rates!(state.heating, state.fluxes, state.atmosphere;
                       gravity = 9.80665, heat_capacity = 1004.0)
    end
    radiation_update!() = begin
        gas!()
        lw!()
        sw!()
        heat!()
    end

    radiation_update!()
    gas_ms, gas_p90 = median_runtime_ms!(gas!, samples)
    lw_ms, lw_p90 = median_runtime_ms!(lw!, samples)
    sw_ms, sw_p90 = median_runtime_ms!(sw!, samples)
    heating_ms, heating_p90 = median_runtime_ms!(heat!, samples)
    total_ms, total_p90 = median_runtime_ms!(radiation_update!, samples)

    return (
        case = "rcemip_style_column_batch",
        date = string(Dates.now()),
        benchmark_status = "scaffold_no_rrtmgp_baseline",
        workload_class = "RCEMIP-style non-spinup column-batch radiation update",
        backend = "CPU",
        device = "host",
        precision = "Float64",
        nx = nx,
        ny = ny,
        columns = length(states),
        layers = nlayers,
        lw_gpoints = 4,
        sw_gpoints = 3,
        simulated_dynamics_timestep_seconds = 300,
        radiation_update_cadence_seconds = 3600,
        simulated_radiation_updates = 1,
        includes_breeze_timestep_overhead = false,
        rrtmgp_baseline_runtime_ms_median = nothing,
        speedup_vs_rrtmgp = nothing,
        nsight_systems_report = nothing,
        nsight_compute_report = nothing,
        gas_optics_runtime_ms_median = gas_ms,
        gas_optics_runtime_ms_p90 = gas_p90,
        longwave_solver_runtime_ms_median = lw_ms,
        longwave_solver_runtime_ms_p90 = lw_p90,
        shortwave_solver_runtime_ms_median = sw_ms,
        shortwave_solver_runtime_ms_p90 = sw_p90,
        heating_runtime_ms_median = heating_ms,
        heating_runtime_ms_p90 = heating_p90,
        radiation_update_runtime_ms_median = total_ms,
        radiation_update_runtime_ms_p90 = total_p90,
        gas_optics_allocations = (@allocated gas!()),
        longwave_solver_allocations = (@allocated lw!()),
        shortwave_solver_allocations = (@allocated sw!()),
        heating_allocations = (@allocated heat!()),
        radiation_update_allocations = (@allocated radiation_update!()),
    )
end

function json_value(value)
    if value === nothing
        return "null"
    elseif value isa AbstractString
        return "\"$value\""
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa NamedTuple
        return json_object(value)
    elseif value isa AbstractVector
        return "[" * join(json_value.(value), ", ") * "]"
    else
        return string(value)
    end
end

function json_object(result)
    fields = propertynames(result)
    lines = String["{"]
    for (i, name) in enumerate(fields)
        value = getproperty(result, name)
        comma = i == length(fields) ? "" : ","
        push!(lines, "  \"$(name)\": $(json_value(value))$(comma)")
    end
    push!(lines, "}")
    return join(lines, "\n")
end

function json_string(results)
    body = join(json_object.(results), ",\n")
    return "{\n  \"cases\": [\n" * replace(body, "\n" => "\n    ") * "\n  ]\n}"
end

function markdown_analytic(result)
    return """
    # Analytic Column Benchmark

    | Metric | Value |
    |---|---:|
    | Case | $(result.case) |
    | Backend | $(result.backend) |
    | Precision | $(result.precision) |
    | Columns | $(result.columns) |
    | Layers | $(result.layers) |
    | Runtime minimum | $(@sprintf("%.6f", result.runtime_ms_minimum)) ms |
    | Runtime median | $(@sprintf("%.6f", result.runtime_ms_median)) ms |
    | Runtime p90 | $(@sprintf("%.6f", result.runtime_ms_p90)) ms |
    | Allocations | $(result.allocations) bytes |
    """
end

function markdown_staged(result)
    return """

    # Staged ecCKD Cloudless Column Benchmark

    | Metric | Value |
    |---|---:|
    | Case | $(result.case) |
    | Backend | $(result.backend) |
    | Precision | $(result.precision) |
    | Columns | $(result.columns) |
    | Layers | $(result.layers) |
    | LW g-points | $(result.lw_gpoints) |
    | SW g-points | $(result.sw_gpoints) |
    | Gas optics median | $(@sprintf("%.6f", result.gas_optics_runtime_ms_median)) ms |
    | Longwave solver median | $(@sprintf("%.6f", result.longwave_solver_runtime_ms_median)) ms |
    | Shortwave solver median | $(@sprintf("%.6f", result.shortwave_solver_runtime_ms_median)) ms |
    | Heating median | $(@sprintf("%.6f", result.heating_runtime_ms_median)) ms |
    | Total median | $(@sprintf("%.6f", result.total_runtime_ms_median)) ms |
    | Gas optics allocations | $(result.gas_optics_allocations) bytes |
    | Longwave solver allocations | $(result.longwave_solver_allocations) bytes |
    | Shortwave solver allocations | $(result.shortwave_solver_allocations) bytes |
    | Heating allocations | $(result.heating_allocations) bytes |
    | Total allocations | $(result.total_allocations) bytes |
    """
end

function markdown_rcemip(result)
    baseline = result.rrtmgp_baseline_runtime_ms_median === nothing ?
        "not yet measured" : string(result.rrtmgp_baseline_runtime_ms_median)
    speedup = result.speedup_vs_rrtmgp === nothing ?
        "not yet measured" : string(result.speedup_vs_rrtmgp)
    return """

    # RCEMIP-Style Column-Batch Benchmark

    This is a local scaffold for the future Breeze+RRTMGP production
    comparison. It runs a non-spinup, nontrivial 3D column-batch radiation
    update with RCEMIP-like moist thermodynamic structure, but it does not yet
    include Breeze timestep overhead, an RRTMGP baseline, GPU execution, or
    Nsight profiles.

    | Metric | Value |
    |---|---:|
    | Case | $(result.case) |
    | Status | $(result.benchmark_status) |
    | Backend | $(result.backend) |
    | Grid | $(result.nx) x $(result.ny) |
    | Columns | $(result.columns) |
    | Layers | $(result.layers) |
    | Radiation cadence | $(result.radiation_update_cadence_seconds) s |
    | Dynamics timestep | $(result.simulated_dynamics_timestep_seconds) s |
    | Gas optics median | $(@sprintf("%.6f", result.gas_optics_runtime_ms_median)) ms |
    | Longwave solver median | $(@sprintf("%.6f", result.longwave_solver_runtime_ms_median)) ms |
    | Shortwave solver median | $(@sprintf("%.6f", result.shortwave_solver_runtime_ms_median)) ms |
    | Heating median | $(@sprintf("%.6f", result.heating_runtime_ms_median)) ms |
    | Radiation update median | $(@sprintf("%.6f", result.radiation_update_runtime_ms_median)) ms |
    | Radiation update allocations | $(result.radiation_update_allocations) bytes |
    | RRTMGP baseline | $(baseline) |
    | Speedup vs RRTMGP | $(speedup) |
    """
end

markdown_string(results) =
    markdown_analytic(results[1]) * "\n" *
    markdown_staged(results[2]) * "\n" *
    markdown_rcemip(results[3])

function main()
    results = [run_case(), run_staged_ecckd_case(), run_rcemip_style_case()]
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    write(joinpath(results_dir, "latest.json"), json_string(results))
    write(joinpath(results_dir, "latest.md"), markdown_string(results))

    println("Analytic column benchmark")
    println("median runtime: $(@sprintf("%.6f", results[1].runtime_ms_median)) ms")
    println("p90 runtime: $(@sprintf("%.6f", results[1].runtime_ms_p90)) ms")
    println("allocations: $(results[1].allocations) bytes")
    println("Staged ecCKD cloudless column benchmark")
    println("gas optics median runtime: $(@sprintf("%.6f", results[2].gas_optics_runtime_ms_median)) ms")
    println("total median runtime: $(@sprintf("%.6f", results[2].total_runtime_ms_median)) ms")
    println("total allocations: $(results[2].total_allocations) bytes")
    println("RCEMIP-style column-batch benchmark")
    println("radiation update median runtime: $(@sprintf("%.6f", results[3].radiation_update_runtime_ms_median)) ms")
    println("radiation update allocations: $(results[3].radiation_update_allocations) bytes")
end

main()
