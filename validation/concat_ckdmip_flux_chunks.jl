using Dates
using NCDatasets

function parse_args(args)
    length(args) >= 2 || error("usage: julia concat_ckdmip_flux_chunks.jl RAW_1.h5 [RAW_2.h5 ...] output.h5")
    return args[1:end-1], args[end]
end

function copy_attributes!(dest, src)
    for key in keys(src.attrib)
        dest.attrib[key] = src.attrib[key]
    end
end

function dim_length(dataset, name)
    return dataset.dim[name]
end

function concatenate_variable!(out, sources, name)
    source_var = sources[1][name]
    dim_names = collect(dimnames(source_var))
    attribs = Dict(key => source_var.attrib[key] for key in keys(source_var.attrib))
    out_var = defVar(out, name, eltype(source_var.var), Tuple(dim_names); attrib = attribs)
    full_selectors = ntuple(_ -> Colon(), length(dim_names))
    column_axis = findfirst(==("column"), dim_names)
    if column_axis === nothing
        out_var[full_selectors...] = source_var[full_selectors...]
        return
    end

    start = 1
    for source in sources
        data = source[name][full_selectors...]
        count = size(data, column_axis)
        selectors = ntuple(axis -> begin
            axis == column_axis ? (start:start + count - 1) : Colon()
        end, ndims(data))
        out_var[selectors...] = data
        start += count
    end
end

function concat_ckdmip_flux_chunks(input_paths, output_path)
    isempty(input_paths) && error("no input chunks supplied")
    sources = NCDataset.(input_paths)
    try
        first_source = sources[1]
        mkpath(dirname(output_path))
        isfile(output_path) && rm(output_path)
        out = NCDataset(output_path, "c")
        try
            for (name, length) in first_source.dim
                defDim(out, name, name == "column" ? Inf : length)
            end
            copy_attributes!(out, first_source)
            previous_history = haskey(first_source.attrib, "history") ? string(first_source.attrib["history"]) : ""
            out.attrib["history"] =
                "$(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")): concat_ckdmip_flux_chunks.jl " *
                join(input_paths, " ") * " " * output_path *
                (isempty(previous_history) ? "" : "\n" * previous_history)
            out.attrib["NCO"] = "not used; concatenated by validation/concat_ckdmip_flux_chunks.jl"

            for name in keys(first_source)
                concatenate_variable!(out, sources, name)
            end
        finally
            close(out)
        end
    finally
        foreach(close, sources)
    end
    return output_path
end

function main(args = ARGS)
    input_paths, output_path = parse_args(args)
    concat_ckdmip_flux_chunks(input_paths, output_path)
    println("Wrote $(output_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
