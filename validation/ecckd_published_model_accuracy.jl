using Dates
using Printf

using NCDatasets

include(joinpath(@__DIR__, "reduced_ecckd_accuracy.jl"))

const PUBLISHED_MODEL_ACCURACY_JSON =
    joinpath(@__DIR__, "results", "ecckd_published_model_accuracy.json")
const PUBLISHED_MODEL_ACCURACY_MD =
    joinpath(@__DIR__, "results", "ecckd_published_model_accuracy.md")

const PUBLISHED_MODEL_SPECS = (
    (
        label = "official ecCKD 1.0 32-LW x 32-SW climate model",
        longwave = :longwave_32,
        shortwave = :shortwave_32,
    ),
    (
        label = "official ecCKD 1.2 64-LW x 64-SW climate model",
        longwave = :longwave_64,
        shortwave = :shortwave_64,
    ),
    (
        label = "official ecCKD 1.2/1.4 64-LW x 96-SW climate/vfine model",
        longwave = :longwave_64,
        shortwave = :shortwave_96,
    ),
)

const PUBLISHED_MODEL_ISOLATION_SPECS = (
    (
        label = "component isolation: published 32-LW x 64-SW",
        longwave = :longwave_32,
        shortwave = :shortwave_64,
        diagnostic = "32-LW reference with 64-SW published component",
    ),
    (
        label = "component isolation: published 32-LW x 96-SW",
        longwave = :longwave_32,
        shortwave = :shortwave_96,
        diagnostic = "32-LW reference with 96-SW published component",
    ),
    (
        label = "component isolation: published 64-LW x 32-SW",
        longwave = :longwave_64,
        shortwave = :shortwave_32,
        diagnostic = "64-LW published component with 32-SW reference",
    ),
)

function published_model(spec)
    return read_ecckd_tabulated_gas_optics(
        official_ecckd_definition_path(spec.longwave),
        official_ecckd_definition_path(spec.shortwave);
        gas_names = OFFICIAL_ECCKD_GASES,
        h2o_mole_fraction = env_float("RH_ECCKD_H2O_MOLE_FRACTION", 0.005),
    )
end

function hard_objective(cases)
    rows = NamedTuple[]
    for case in cases
        for variable in (:lw_up, :lw_down, :sw_up, :sw_down)
            metrics = getproperty(case.variables, variable)
            push!(rows, (
                case = case.case,
                metric = "$(variable)_rmse",
                value = metrics.rmse,
                threshold = ACCEPTANCE_THRESHOLDS.flux_rmse_w_m2,
                normalized_value = metrics.rmse / ACCEPTANCE_THRESHOLDS.flux_rmse_w_m2,
            ))
            push!(rows, (
                case = case.case,
                metric = "$(variable)_max_abs",
                value = metrics.max_abs,
                threshold = ACCEPTANCE_THRESHOLDS.flux_max_abs_w_m2,
                normalized_value = metrics.max_abs / ACCEPTANCE_THRESHOLDS.flux_max_abs_w_m2,
            ))
        end
        heating = case.variables.heating_rate
        push!(rows, (
            case = case.case,
            metric = "heating_rate_rmse",
            value = heating.rmse,
            threshold = ACCEPTANCE_THRESHOLDS.heating_rate_rmse_k_day,
            normalized_value = heating.rmse / ACCEPTANCE_THRESHOLDS.heating_rate_rmse_k_day,
        ))
        push!(rows, (
            case = case.case,
            metric = "heating_rate_max_abs",
            value = heating.max_abs,
            threshold = ACCEPTANCE_THRESHOLDS.heating_rate_max_abs_k_day,
            normalized_value = heating.max_abs / ACCEPTANCE_THRESHOLDS.heating_rate_max_abs_k_day,
        ))
        push!(rows, (
            case = case.case,
            metric = "toa_forcing",
            value = case.toa_forcing_max_abs,
            threshold = ACCEPTANCE_THRESHOLDS.toa_forcing_abs_error_w_m2,
            normalized_value =
                case.toa_forcing_max_abs / ACCEPTANCE_THRESHOLDS.toa_forcing_abs_error_w_m2,
        ))
        push!(rows, (
            case = case.case,
            metric = "surface_forcing",
            value = case.surface_forcing_max_abs,
            threshold = ACCEPTANCE_THRESHOLDS.surface_forcing_abs_error_w_m2,
            normalized_value =
                case.surface_forcing_max_abs /
                ACCEPTANCE_THRESHOLDS.surface_forcing_abs_error_w_m2,
        ))
    end
    worst = argmax(row -> row.normalized_value, rows)
    return (
        value = worst.normalized_value,
        case = worst.case,
        metric = worst.metric,
        metric_value = worst.value,
        threshold = worst.threshold,
    )
end

function published_model_result(spec)
    model = published_model(spec)
    cases = [case_metrics(case, model) for case in REDUCED_CASES]
    objective = hard_objective(cases)
    return (
        label = spec.label,
        diagnostic = hasproperty(spec, :diagnostic) ? spec.diagnostic : "",
        longwave_key = String(spec.longwave),
        shortwave_key = String(spec.shortwave),
        ng_lw = size(model.longwave_absorption, 1),
        ng_sw = size(model.shortwave_absorption, 1),
        total_gpoints = size(model.longwave_absorption, 1) +
                        size(model.shortwave_absorption, 1),
        longwave_path = official_ecckd_definition_path(spec.longwave),
        shortwave_path = official_ecckd_definition_path(spec.shortwave),
        passed_hard_thresholds = all(case -> case.passed_hard_thresholds, cases),
        worst_toa_forcing_abs_error_w_m2 =
            maximum(case -> case.toa_forcing_max_abs, cases),
        worst_surface_forcing_abs_error_w_m2 =
            maximum(case -> case.surface_forcing_max_abs, cases),
        hard_objective = objective,
        cases = cases,
    )
end

function published_model_accuracy()
    ENV["RH_CANDIDATE_GAS_OPTICS"] = "official_ecckd"
    models = [published_model_result(spec) for spec in PUBLISHED_MODEL_SPECS]
    isolation_diagnostics =
        [published_model_result(spec) for spec in PUBLISHED_MODEL_ISOLATION_SPECS]
    return (
        case = "ecckd_published_model_accuracy",
        timestamp_utc = string(Dates.now()),
        status = all(model -> model.passed_hard_thresholds, models) ? "passed" :
                 "failed_threshold",
        reference_scope = collect(REDUCED_CASE_NAMES),
        acceptance_thresholds = ACCEPTANCE_THRESHOLDS,
        models = models,
        isolation_diagnostics = isolation_diagnostics,
    )
end

function markdown_report(result)
    lines = String[
        "# Published ecCKD Model Accuracy",
        "",
        "Status: **$(result.status)**",
        "",
        "Reference scope: clean ecCKD cloudless/no-aerosol tropical and RCEMIP-style cases.",
        "",
        "| Model | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for model in result.models
        push!(lines,
              "| $(model.label) | $(model.ng_lw) | $(model.ng_sw) | $(model.passed_hard_thresholds) | $(@sprintf("%.12g", model.worst_toa_forcing_abs_error_w_m2)) W m^-2 | $(@sprintf("%.12g", model.worst_surface_forcing_abs_error_w_m2)) W m^-2 | $(@sprintf("%.12g", model.hard_objective.value)) | $(model.hard_objective.metric) |")
    end
    append!(lines, [
        "",
        "## Mixed-component isolation diagnostics",
        "",
        "| Diagnostic | LW | SW | Passed | Worst TOA forcing | Worst surface forcing | Hard objective | Limiting metric |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ])
    for model in result.isolation_diagnostics
        push!(lines,
              "| $(model.diagnostic) | $(model.ng_lw) | $(model.ng_sw) | $(model.passed_hard_thresholds) | $(@sprintf("%.12g", model.worst_toa_forcing_abs_error_w_m2)) W m^-2 | $(@sprintf("%.12g", model.worst_surface_forcing_abs_error_w_m2)) W m^-2 | $(@sprintf("%.12g", model.hard_objective.value)) | $(model.hard_objective.metric) |")
    end
    append!(lines, [
        "",
        "This artifact evaluates published full-accuracy ecCKD CKD-definition combinations directly against the same clean package-native ecRad reference cases used by the reduced-model gate.",
        "The mixed-component rows are diagnostic only: they are not published ecCKD products, but they localize the current 64/96 parity failure by swapping one published component at a time against the passing 32x32 baseline.",
    ])
    return join(lines, "\n") * "\n"
end

function main()
    result = published_model_accuracy()
    mkpath(dirname(PUBLISHED_MODEL_ACCURACY_JSON))
    write(PUBLISHED_MODEL_ACCURACY_JSON, json_object(result) * "\n")
    write(PUBLISHED_MODEL_ACCURACY_MD, markdown_report(result))
    print(markdown_report(result))
    println("Wrote $PUBLISHED_MODEL_ACCURACY_JSON")
    println("Wrote $PUBLISHED_MODEL_ACCURACY_MD")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
