using Dates

const ABR_ROOT = normpath(joinpath(@__DIR__, ".."))
if Base.find_package("Lightflux") === nothing
    push!(LOAD_PATH, ABR_ROOT)
end

using Lightflux

const TRAINING_MANIFEST_JSON =
    joinpath(@__DIR__, "results", "ecckd_published_training_manifest.json")
const TRAINING_MANIFEST_MD =
    joinpath(@__DIR__, "results", "ecckd_published_training_manifest.md")

function json_escape(text)
    return replace(text, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
end

function json_value(value)
    if value isa AbstractString
        return "\"" * json_escape(value) * "\""
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa NamedTuple
        return json_object(value)
    elseif value isa AbstractVector || value isa Tuple
        return "[" * join(json_value.(value), ", ") * "]"
    elseif value === nothing
        return "null"
    else
        return string(value)
    end
end

function json_object(object)
    names = propertynames(object)
    lines = ["{"]
    for (i, name) in enumerate(names)
        comma = i == length(names) ? "" : ","
        push!(lines, "  \"$(name)\": $(json_value(getproperty(object, name)))$(comma)")
    end
    push!(lines, "}")
    return join(lines, "\n")
end

function active_assignment_lines(path, variable)
    lines = String[]
    for line in split(read(path, String), '\n')
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        occursin(variable * "=", stripped) && push!(lines, stripped)
    end
    return lines
end

function config_value(path, variable)
    prefix = variable * "="
    for line in split(read(path, String), '\n')
        stripped = strip(line)
        startswith(stripped, prefix) || continue
        return strip(chop(stripped; head = length(prefix), tail = 0))
    end
    return nothing
end

function parse_mode_blocks(path)
    blocks = NamedTuple[]
    current_names = String[]
    current_lines = String[]
    in_block = false
    for line in split(read(path, String), '\n')
        stripped = strip(line)
        m = match(r"^([A-Za-z0-9_-]+(?:\s*\|\s*[A-Za-z0-9_-]+)*)\)$", stripped)
        if m !== nothing
            current_names = [strip(name) for name in split(m.captures[1], "|")]
            current_lines = String[]
            in_block = true
            continue
        end
        if in_block && stripped == ";;"
            assignments = logical_assignments(current_lines)
            push!(blocks, (
                names = current_names,
                assignments = assignments,
            ))
            current_names = String[]
            current_lines = String[]
            in_block = false
            continue
        end
        in_block && push!(current_lines, line)
    end
    return blocks
end

function quote_count(text)
    return count(==('"'), text)
end

function logical_assignments(lines)
    variables = ("TRAINING", "GASLIST", "SPECIFIC_OPTIONS", "EXTRA_ARGS")
    assignments = String[]
    buffer = nothing
    skipping_training_both_branch = false
    for line in lines
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        if skipping_training_both_branch
            stripped == "fi" && (skipping_training_both_branch = false)
            continue
        end
        if startswith(stripped, "if ") && occursin("TRAINING_BOTH", stripped)
            skipping_training_both_branch = true
            continue
        end
        if buffer !== nothing
            buffer *= " " * stripped
            if iseven(quote_count(buffer))
                push!(assignments, buffer)
                buffer = nothing
            end
            continue
        end
        any(var -> startswith(stripped, var * "="), variables) || continue
        if isodd(quote_count(stripped))
            buffer = stripped
        else
            push!(assignments, stripped)
        end
    end
    buffer === nothing || push!(assignments, buffer)
    return assignments
end

function script_summary(root, relative_path)
    path = joinpath(root, relative_path)
    common_options = active_assignment_lines(path, "COMMON_OPTIONS")
    return (
        path = relative_path,
        exists = isfile(path),
        active_common_options = common_options,
        selected_common_options = isempty(common_options) ? nothing : last(common_options),
        modes = parse_mode_blocks(path),
    )
end

function run_ecckd_published_training_manifest()
    root = ecckd_source_path(require = true)
    config_path = joinpath(root, "test", "config.h")
    required_files = [
        "src/ecckd/optimize_lut.cpp",
        "src/ecckd/solve_adept.cpp",
        "src/ecckd/calc_cost_function_lw.cpp",
        "src/ecckd/calc_cost_function_sw.cpp",
        "test/config.h",
        "test/find_g_points_lw.sh",
        "test/find_g_points_sw.sh",
        "test/create_lut_lw.sh",
        "test/create_lut_sw.sh",
        "test/optimize_lut_lw.sh",
        "test/optimize_lut_sw.sh",
        "test/run_ckd_lw.sh",
        "test/run_ckd_sw.sh",
    ]
    source_files = [(path = file, exists = isfile(joinpath(root, file))) for file in required_files]
    scripts = [
        script_summary(root, "test/optimize_lut_lw.sh"),
        script_summary(root, "test/optimize_lut_sw.sh"),
    ]
    config = (
        ckdmip_data_dir = config_value(config_path, "CKDMIP_DATA_DIR"),
        training_code = config_value(config_path, "TRAINING_CODE"),
        evaluation_code = config_value(config_path, "EVALUATION_CODE"),
        training_both = config_value(config_path, "TRAINING_BOTH"),
        mmm_code = config_value(config_path, "MMM_CODE"),
        idealized_code = config_value(config_path, "IDEALIZED_CODE"),
    )
    passed = all(file -> file.exists, source_files) &&
             all(script -> script.exists && script.selected_common_options !== nothing, scripts)
    return (
        case = "ecckd_published_training_manifest",
        timestamp_utc = string(Dates.now()),
        status = passed ? "passed" : "failed",
        ecckd_source_root = root,
        objective_contract =
            "Keep these ecCKD source files, script modes, data selections, gases, weights, regularization, and stopping settings fixed; replace only the optimizer implementation/settings in the Julia Reactant/Enzyme recovery pipeline.",
        config = config,
        source_files = source_files,
        optimization_scripts = scripts,
        remaining_external_data_requirement =
            "Full CKDMIP line-by-line spectral absorption and LBL flux database under RH_CKDMIP_DATA_PATH.",
    )
end

function markdown_manifest(result)
    lines = String[
        "# ecCKD Published Training Manifest",
        "",
        "Status: **$(result.status)**",
        "",
        "ecCKD source root: `$(result.ecckd_source_root)`",
        "",
        result.objective_contract,
        "",
        "## Config",
        "",
        "| Field | Value |",
        "|---|---|",
        "| CKDMIP data dir template | `$(result.config.ckdmip_data_dir)` |",
        "| Training dataset | `$(result.config.training_code)` |",
        "| Evaluation dataset | `$(result.config.evaluation_code)` |",
        "| Train with both evaluation sets | `$(result.config.training_both)` |",
        "| MMM dataset | `$(result.config.mmm_code)` |",
        "| Idealized dataset | `$(result.config.idealized_code)` |",
        "",
        "## Source Files",
        "",
        "| Path | Present |",
        "|---|---:|",
    ]
    for file in result.source_files
        push!(lines, "| `$(file.path)` | $(file.exists) |")
    end
    push!(lines, "", "## Optimization Scripts", "")
    for script in result.optimization_scripts
        push!(lines, "### `$(script.path)`", "")
        push!(lines, "- Selected common options: `$(script.selected_common_options)`")
        mode_names = [join(mode.names, " | ") for mode in script.modes]
        push!(lines, "- Modes: $(join(mode_names, ", "))", "")
    end
    push!(lines, "Remaining external data requirement: $(result.remaining_external_data_requirement)")
    return join(lines, "\n") * "\n"
end

function ecckd_published_training_manifest_main()
    result = run_ecckd_published_training_manifest()
    mkpath(dirname(TRAINING_MANIFEST_JSON))
    write(TRAINING_MANIFEST_JSON, json_object(result) * "\n")
    write(TRAINING_MANIFEST_MD, markdown_manifest(result))
    print(markdown_manifest(result))
    println("Wrote $TRAINING_MANIFEST_JSON")
    println("Wrote $TRAINING_MANIFEST_MD")
end

if abspath(PROGRAM_FILE) == @__FILE__
    ecckd_published_training_manifest_main()
end
