using Dates
using Printf

include(joinpath(@__DIR__, "reduced_ecckd_optimization_preflight.jl"))

const OFFICIAL_TRAINING_JSON = joinpath(@__DIR__, "results", "official_ecckd_training.json")
const OFFICIAL_TRAINING_MD = joinpath(@__DIR__, "results", "official_ecckd_training.md")

function official_training_status(preflight)
    improved = preflight.directional_search.final_objective <
               preflight.directional_search.initial_objective
    checks_passed = preflight.enzyme_check.status == "passed" &&
                    preflight.reactant_check.status == "passed"
    return improved && checks_passed ? "passed" : "failed"
end

function official_training_report(preflight)
    status = official_training_status(preflight)
    initial_objective = preflight.directional_search.initial_objective
    final_objective = preflight.directional_search.final_objective
    objective_target = preflight.objective_target
    objective_ratio = final_objective / objective_target
    return (
        case = "official_reduced_ecckd_gas_optics_training",
        timestamp_utc = string(Dates.now()),
        status = status,
        reference_scope = preflight.reference_scope,
        training_source = "official ecCKD tabulated gas-optics reduced hard-gate objective",
        parameterization = preflight.parameterization,
        parameter_count = preflight.parameter_count,
        ng_sw = preflight.ng_sw,
        optimizer = preflight.directional_search.method,
        iterations = preflight.directional_search.iterations,
        initial_objective = initial_objective,
        final_objective = final_objective,
        objective_reduction = initial_objective - final_objective,
        objective_ratio = final_objective / initial_objective,
        objective_target = objective_target,
        final_objective_target_ratio = objective_ratio,
        hard_accuracy_target_met = final_objective <= objective_target,
        acceptance_gap_status = preflight.acceptance_gap_status,
        loss_history = [
            step.objective for step in preflight.directional_search.steps
        ],
        reactant_check = preflight.reactant_check,
        enzyme_check = preflight.enzyme_check,
        topology_candidate_scan = preflight.topology_candidate_scan,
        next_required_work = preflight.next_required_work,
        notes = "This is the official/reduced ecCKD training-path artifact: it demonstrates objective construction, trainable parameters, Reactant/Enzyme checks, and deterministic objective reduction on official ecCKD references. It does not close the reduced-accuracy gate until final_objective_target_ratio <= 1.",
    )
end

function markdown_report(result)
    lines = String[
        "# Official/Reduced ecCKD Gas-Optics Training",
        "",
        "Status: **$(result.status)**",
        "",
        "| Field | Value |",
        "|---|---:|",
        "| Trainable shortwave g-points | $(result.ng_sw) |",
        "| Parameter count | $(result.parameter_count) |",
        "| Iterations | $(result.iterations) |",
        "| Initial objective | $(@sprintf("%.12g", result.initial_objective)) |",
        "| Final objective | $(@sprintf("%.12g", result.final_objective)) |",
        "| Objective reduction | $(@sprintf("%.12g", result.objective_reduction)) |",
        "| Final objective / initial objective | $(@sprintf("%.12g", result.objective_ratio)) |",
        "| Final objective / hard target | $(@sprintf("%.12g", result.final_objective_target_ratio)) |",
        "| Hard accuracy target met | $(result.hard_accuracy_target_met) |",
        "| Reactant check | $(result.reactant_check.status) |",
        "| Enzyme check | $(result.enzyme_check.status) |",
        "",
        result.notes,
        "",
        "Next required work: $(result.next_required_work)",
    ]
    return join(lines, "\n") * "\n"
end

function main()
    preflight = reduced_optimization_preflight()
    result = official_training_report(preflight)
    mkpath(dirname(OFFICIAL_TRAINING_JSON))
    write(OFFICIAL_TRAINING_JSON, json_object(result) * "\n")
    write(OFFICIAL_TRAINING_MD, markdown_report(result))
    print(markdown_report(result))
    println("Wrote $OFFICIAL_TRAINING_JSON")
    println("Wrote $OFFICIAL_TRAINING_MD")
    result.status == "passed" || error("official/reduced ecCKD training failed")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
