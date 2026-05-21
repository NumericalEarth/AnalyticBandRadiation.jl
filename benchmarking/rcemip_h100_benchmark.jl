"""
RCEMIP-style H100 benchmark for RadiativeHeating.jl (via the BreezeRadiative-
HeatingExt streaming gas-optics + transport kernel) against the RRTMGP.jl
baseline (via the BreezeRRTMGPExt clear-sky kernel), measured under Breeze's
`update_radiation!` call surface.

Configured by environment variables (defaults shown):

  RADIATIVE_HEATING_BACKEND               "CPU" | "H100"   ("CPU")
  RADIATIVE_HEATING_GAS_MODEL_SOURCE      "validated_ecCKD" |
                                          "synthetic_fixed_coefficients"
                                                            ("validated_ecCKD")
  RADIATIVE_HEATING_ECCKD_LW_PATH         absolute path to the ecCKD LW
                                          ckd-definition.nc to load when
                                          gas-model source is validated_ecCKD
  RADIATIVE_HEATING_ECCKD_SW_PATH         absolute path to the ecCKD SW
                                          ckd-definition.nc to load when
                                          gas-model source is validated_ecCKD
  RADIATIVE_HEATING_NG_LW                 g-points (synthetic path)        (32)
  RADIATIVE_HEATING_NG_SW                 g-points (synthetic path)        (16)
  RADIATIVE_HEATING_NX, NY, NZ            grid dimensions       (32, 32, 64)
  RADIATIVE_HEATING_SAMPLES               benchmark samples (post-warmup)   (5)
  RADIATIVE_HEATING_RUN_RRTMGP            "true" | "false"               ("true")
  RADIATIVE_HEATING_OUTPUT_DIR            output directory   (benchmarking/results)
  RADIATIVE_HEATING_LABEL                 label suffix on the output dir    ("")

Writes `radiative_heating_rcemip_latest.{json,md}` under
`RADIATIVE_HEATING_OUTPUT_DIR`. Schema matches the dedicated Breeze checkout's
benchmark, so `figures/make_pr_figures.jl` can pick the result up directly.
"""

using AnalyticBandRadiation
using Breeze
using CUDA
using ClimaComms
using Dates
using JSON
using NCDatasets
using Oceananigans
using Printf
using RRTMGP
using Statistics

using Breeze.AtmosphereModels: update_radiation!, RadiativeHeatingOptics,
                                RadiativeTransferModel, ClearSkyOptics
using Breeze.Thermodynamics: ThermodynamicConstants
using Oceananigans.Architectures: CPU, GPU

const DAY = 24 * 60 * 60
const OFFICIAL_ECCKD_GASES = (:composite, :h2o, :o3, :co2, :ch4, :n2o, :cfc11, :cfc12)

envint(name, default) = parse(Int, get(ENV, name, string(default)))
envfloat(name, default) = parse(Float64, get(ENV, name, string(default)))
envbool(name, default) = lowercase(get(ENV, name, string(default))) in ("true", "1", "yes", "on")

function backend_name()
    requested = uppercase(get(ENV, "RADIATIVE_HEATING_BACKEND", "CPU"))
    if requested == "H100"
        CUDA.functional() || error("RADIATIVE_HEATING_BACKEND=H100 but CUDA is not functional")
        return "H100", GPU()
    end
    return "CPU", CPU()
end

function rcemip_grid(arch, FT, Nx, Ny, Nz)
    return RectilinearGrid(arch, FT;
                           size = (Nx, Ny, Nz),
                           x = (0, 96e3),
                           y = (-12, 12),
                           z = (0, 30e3),
                           topology = (Periodic, Periodic, Bounded))
end

function make_synthetic_gas_optics(FT; ng_lw = 32, ng_sw = 16)
    lw = zeros(FT, ng_lw, 2)
    sw = zeros(FT, ng_sw, 2)
    for ig in 1:ng_lw
        lw[ig, 1] = FT(16 + 0.25ig)
        lw[ig, 2] = FT(18 + 0.1ig)
    end
    for ig in 1:ng_sw
        sw[ig, 1] = FT(220 + 0.8ig)
        sw[ig, 2] = FT(6 + 0.2ig)
    end
    return EcCKDGasOpticsModel(
        gas_names = (:h2o, :co2),
        longwave_absorption = lw,
        shortwave_absorption = sw,
        longwave_weights = fill(inv(FT(ng_lw)), ng_lw),
        shortwave_weights = fill(inv(FT(ng_sw)), ng_sw),
        longwave_source_scale = fill(one(FT), ng_lw),
    )
end

function official_ecckd_gas_optics()
    haskey(ENV, "RADIATIVE_HEATING_ECCKD_LW_PATH") ||
        error("RADIATIVE_HEATING_ECCKD_LW_PATH must point at an ecCKD LW definition.nc when GAS_MODEL_SOURCE=validated_ecCKD")
    haskey(ENV, "RADIATIVE_HEATING_ECCKD_SW_PATH") ||
        error("RADIATIVE_HEATING_ECCKD_SW_PATH must point at an ecCKD SW definition.nc when GAS_MODEL_SOURCE=validated_ecCKD")
    lw_path = ENV["RADIATIVE_HEATING_ECCKD_LW_PATH"]
    sw_path = ENV["RADIATIVE_HEATING_ECCKD_SW_PATH"]
    return read_ecckd_tabulated_gas_optics(lw_path, sw_path;
                                           gas_names = OFFICIAL_ECCKD_GASES,
                                           h2o_mole_fraction = envfloat("RADIATIVE_HEATING_ECCKD_H2O_MOLE_FRACTION", 0.005))
end

function benchmark_gas_model(FT)
    source = get(ENV, "RADIATIVE_HEATING_GAS_MODEL_SOURCE", "validated_ecCKD")
    if source == "validated_ecCKD"
        gas_optics = official_ecckd_gas_optics()
        return (
            gas_optics = gas_optics,
            kind = "official_ecCKD_$(size(gas_optics.longwave_absorption, 1))_lw_$(size(gas_optics.shortwave_absorption, 1))_sw",
            source = "validated_ecCKD",
        )
    elseif source == "synthetic_fixed_coefficients"
        ng_lw = envint("RADIATIVE_HEATING_NG_LW", 32)
        ng_sw = envint("RADIATIVE_HEATING_NG_SW", 16)
        return (
            gas_optics = make_synthetic_gas_optics(FT; ng_lw, ng_sw),
            kind = "fixed_ecCKD_$(ng_lw)_lw_$(ng_sw)_sw",
            source = "synthetic_fixed_coefficients",
        )
    end
    error("unsupported RADIATIVE_HEATING_GAS_MODEL_SOURCE=$source")
end

function default_gas_values(FT, source)
    if source == "validated_ecCKD"
        return Dict{Symbol, FT}(
            :co2 => FT(envfloat("RADIATIVE_HEATING_CO2_VMR", 420.0e-6)),
            :o3 => FT(envfloat("RADIATIVE_HEATING_O3_VMR", 0.0)),
            :ch4 => FT(envfloat("RADIATIVE_HEATING_CH4_VMR", 1.8e-6)),
            :n2o => FT(envfloat("RADIATIVE_HEATING_N2O_VMR", 330.0e-9)),
            :cfc11 => FT(envfloat("RADIATIVE_HEATING_CFC11_VMR", 230.0e-12)),
            :cfc12 => FT(envfloat("RADIATIVE_HEATING_CFC12_VMR", 520.0e-12)),
        )
    end
    return Dict{Symbol, FT}(:co2 => FT(envfloat("RADIATIVE_HEATING_CO2_VMR", 400.0e-6)))
end

function build_model(radiation, grid, constants)
    reference_state = ReferenceState(grid, constants;
                                     surface_pressure = 101325,
                                     potential_temperature = 300)
    dynamics = AnelasticDynamics(reference_state)
    clock = Clock(time = DateTime(2024, 7, 15, 12, 0, 0))
    model = AtmosphereModel(grid; clock, dynamics,
                            formulation = :LiquidIcePotentialTemperature,
                            radiation)
    θ(x, y, z) = 300 + 18 * (z / 30e3)^1.25
    qᵗ(x, y, z) = 0.018 * exp(-z / 2400)
    set!(model; θ, qᵗ)
    return model
end

function median_seconds!(f, samples)
    times = Float64[]
    for _ in 1:samples
        GC.gc(false)
        push!(times, @elapsed f())
    end
    return median(times)
end

sync_backend(backend) = backend == "H100" ? CUDA.synchronize() : nothing

function main()
    FT = Float64
    Nx = envint("RADIATIVE_HEATING_NX", 32)
    Ny = envint("RADIATIVE_HEATING_NY", 32)
    Nz = envint("RADIATIVE_HEATING_NZ", 64)
    samples = envint("RADIATIVE_HEATING_SAMPLES", 5)
    run_rrtmgp = envbool("RADIATIVE_HEATING_RUN_RRTMGP", true)
    label_suffix = get(ENV, "RADIATIVE_HEATING_LABEL", "")
    output_dir = get(ENV, "RADIATIVE_HEATING_OUTPUT_DIR",
                     joinpath(@__DIR__, "results",
                              "rcemip_$(Nx)x$(Ny)x$(Nz)" *
                              (isempty(label_suffix) ? "" : "_$(label_suffix)")))
    mkpath(output_dir)

    backend, arch = backend_name()
    grid = rcemip_grid(arch, FT, Nx, Ny, Nz)
    constants = Breeze.Thermodynamics.ThermodynamicConstants()

    @info "Benchmark setup" backend Nx Ny Nz samples columns=Nx*Ny output_dir

    gas_model = benchmark_gas_model(FT)
    gas_values = default_gas_values(FT, gas_model.source)

    # ---- RadiativeHeating path ----
    rh_seconds = Inf
    rh_supported = true
    rh_error = ""
    try
        rh_radiation = Base.invokelatest(() ->
            RadiativeTransferModel(grid, RadiativeHeatingOptics(), constants;
                                   gas_optics = gas_model.gas_optics,
                                   gas_values,
                                   surface_temperature = 300,
                                   surface_albedo = 0.07,
                                   surface_emissivity = 0.98,
                                   solar_constant = 551))
        rh_model = build_model(rh_radiation, grid, constants)
        Base.invokelatest(update_radiation!, rh_model.radiation, rh_model)
        sync_backend(backend)
        @info "RadiativeHeating warmup complete; timing $(samples) samples"
        rh_seconds = median_seconds!(samples) do
            Base.invokelatest(update_radiation!, rh_model.radiation, rh_model)
            sync_backend(backend)
        end
        @info "RadiativeHeating median update" ms = 1000 * rh_seconds
    catch err
        rh_supported = false
        rh_error = sprint(showerror, err)
        @warn "RadiativeHeating path failed" err
    end

    # ---- RRTMGP baseline ----
    rrtmgp_seconds = Inf
    rrtmgp_supported = true
    rrtmgp_error = ""
    if run_rrtmgp
        try
            rrtmgp_radiation = Base.invokelatest(() ->
                RadiativeTransferModel(grid, ClearSkyOptics(), constants;
                                       surface_temperature = 300,
                                       surface_emissivity = 0.98,
                                       surface_albedo = 0.07,
                                       solar_constant = 551,
                                       coordinate = (0.0, 0.0)))
            rrtmgp_model = build_model(rrtmgp_radiation, grid, constants)
            Base.invokelatest(update_radiation!, rrtmgp_model.radiation, rrtmgp_model)
            sync_backend(backend)
            @info "RRTMGP warmup complete; timing $(samples) samples"
            rrtmgp_seconds = median_seconds!(samples) do
                Base.invokelatest(update_radiation!, rrtmgp_model.radiation, rrtmgp_model)
                sync_backend(backend)
            end
            @info "RRTMGP median update" ms = 1000 * rrtmgp_seconds
        catch err
            rrtmgp_supported = false
            rrtmgp_error = sprint(showerror, err)
            @warn "RRTMGP path failed" err
        end
    end

    speedup = (isfinite(rh_seconds) && isfinite(rrtmgp_seconds)) ?
              rrtmgp_seconds / rh_seconds : 0.0

    result = Dict(
        "case" => "radiative_heating_rcemip_benchmark",
        "timestamp_utc" => string(now()),
        "backend" => backend,
        "gas_model_kind" => gas_model.kind,
        "gas_model_source" => gas_model.source,
        "ng_lw" => size(gas_model.gas_optics.longwave_absorption, 1),
        "ng_sw" => size(gas_model.gas_optics.shortwave_absorption, 1),
        "grid" => Dict(
            "nx" => Nx, "ny" => Ny, "nz" => Nz,
            "columns" => Nx * Ny,
        ),
        "samples" => samples,
        "radiative_heating_runtime_supported" => rh_supported,
        "radiative_heating_update_median_ms" => isfinite(rh_seconds) ? 1000 * rh_seconds : nothing,
        "radiative_heating_error" => rh_error,
        "rrtmgp_runtime_supported" => rrtmgp_supported,
        "rrtmgp_update_median_ms" => isfinite(rrtmgp_seconds) ? 1000 * rrtmgp_seconds : nothing,
        "rrtmgp_error" => rrtmgp_error,
        "radiation_update_speedup" => speedup,
        "gas_values" => Dict(string(k) => v for (k, v) in gas_values),
    )

    json_path = joinpath(output_dir, "radiative_heating_rcemip_latest.json")
    open(json_path, "w") do io
        JSON.print(io, result, 2)
    end

    md_path = joinpath(output_dir, "radiative_heating_rcemip_latest.md")
    open(md_path, "w") do io
        println(io, "# RadiativeHeating RCEMIP benchmark")
        println(io)
        println(io, "- backend: $(result["backend"])")
        println(io, "- grid: $(Nx) × $(Ny) × $(Nz) ($(Nx * Ny) columns)")
        println(io, "- samples: $samples (post-warmup median)")
        println(io, "- gas model: $(result["gas_model_kind"]) ($(result["gas_model_source"]))")
        println(io, "- ng_lw / ng_sw: $(result["ng_lw"]) / $(result["ng_sw"])")
        println(io)
        if result["radiative_heating_update_median_ms"] !== nothing
            @printf(io, "- RadiativeHeating update median: %.3f ms\n",
                    result["radiative_heating_update_median_ms"])
        else
            println(io, "- RadiativeHeating update: failed ($(rh_error))")
        end
        if result["rrtmgp_update_median_ms"] !== nothing
            @printf(io, "- RRTMGP update median: %.3f ms\n", result["rrtmgp_update_median_ms"])
        elseif run_rrtmgp
            println(io, "- RRTMGP update: failed ($(rrtmgp_error))")
        else
            println(io, "- RRTMGP update: skipped")
        end
        if speedup > 0
            @printf(io, "- speedup vs RRTMGP: %.2fx\n", speedup)
        end
    end

    @info "Wrote benchmark artifact" json = json_path md = md_path

    return result
end

main()
