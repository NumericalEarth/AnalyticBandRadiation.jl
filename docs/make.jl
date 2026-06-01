using Documenter
using NumericalRadiation

makedocs(
    sitename = "NumericalRadiation.jl",
    modules = [NumericalRadiation],
    authors = "NumericalEarth organization and contributors",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://NumericalEarth.github.io/AnalyticBandRadiation.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Architecture" => "architecture.md",
        "Gas optics" => [
            "ecCKD files" => "gas_optics/ecckd_files.md",
            "CKDMIP training data" => "gas_optics/ckdmip_training_data.md",
        ],
        "Longwave" => "longwave.md",
        "Shortwave" => "shortwave.md",
        "Single-column examples" => "single_column.md",
        "Notation" => "notation.md",
        "API reference" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/NumericalEarth/AnalyticBandRadiation.jl.git",
    devbranch = "main",
    push_preview = true,
)
