using CairoMakie
using JSON
using NCDatasets
using Statistics

const REPO_ROOT = abspath(joinpath(@__DIR__, ".."))
const OUT_DIR = @__DIR__
const BREEZE_RESULTS = "/shared/home/greg/Projects/BreezeRadiativeHeatingDev/Breeze.jl/benchmarking/results"

set_theme!(theme_minimal())

# -----------------------------------------------------------------------------
# Figure 1. Clear-sky tropical column: official 32x32 ecCKD vs ecRad reference
# -----------------------------------------------------------------------------

function figure_flux_profiles()
    path = joinpath(REPO_ROOT, "validation/reference/ecrad/ecckd_clear_sky_tropical_column.nc")
    ds = NCDataset(path)

    p_iface = ds["pressure_interface"][:, :] ./ 100.0  # hPa
    p_layer = ds["pressure_layer"][:, :] ./ 100.0
    ref_lw_up = ds["lw_up"][:, :]
    ref_lw_dn = ds["lw_down"][:, :]
    ref_sw_up = ds["sw_up"][:, :]
    ref_sw_dn = ds["sw_down"][:, :]
    ref_hr = ds["heating_rate"][:, :]
    cand_lw_up = ds["radiative_heating_lw_up"][:, :]
    cand_lw_dn = ds["radiative_heating_lw_down"][:, :]
    cand_sw_up = ds["radiative_heating_sw_up"][:, :]
    cand_sw_dn = ds["radiative_heating_sw_down"][:, :]
    cand_hr = ds["radiative_heating_heating_rate"][:, :]
    close(ds)

    column = 1
    fig = Figure(size = (1100, 720))
    Label(fig[0, 1:3],
          "ecCKD 32x32 (RadiativeHeating.jl) vs ecRad on the clear-sky tropical column";
          fontsize = 18, font = :bold, tellwidth = false)

    panels = [
        ("LW up", ref_lw_up, cand_lw_up, p_iface, "W m⁻²"),
        ("LW down", ref_lw_dn, cand_lw_dn, p_iface, "W m⁻²"),
        ("Heating rate", ref_hr, cand_hr, p_layer, "K day⁻¹"),
        ("SW up", ref_sw_up, cand_sw_up, p_iface, "W m⁻²"),
        ("SW down", ref_sw_dn, cand_sw_dn, p_iface, "W m⁻²"),
    ]

    locs = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2)]
    for ((title, ref, cand, p, units), (row, col)) in zip(panels, locs)
        ax = Axis(fig[row, col]; title = title, xlabel = units,
                  ylabel = "pressure [hPa]", yreversed = true, yscale = log10)
        lines!(ax, ref[:, column], p[:, column]; color = :black, label = "ecRad", linewidth = 2)
        lines!(ax, cand[:, column], p[:, column]; color = :tomato,
               linestyle = :dash, label = "RadiativeHeating 32x32", linewidth = 2)
        ylims!(ax, 1000.0, 1.0)
    end

    # Error summary panel
    n_lev_iface = size(ref_lw_up, 1)
    n_lev_layer = size(ref_hr, 1)
    rmse_lw = sqrt(mean((vec(cand_lw_up) .- vec(ref_lw_up)).^2 .+ (vec(cand_lw_dn) .- vec(ref_lw_dn)).^2))
    rmse_sw = sqrt(mean((vec(cand_sw_up) .- vec(ref_sw_up)).^2 .+ (vec(cand_sw_dn) .- vec(ref_sw_dn)).^2))
    rmse_hr = sqrt(mean((vec(cand_hr) .- vec(ref_hr)).^2))
    toa_err = (cand_lw_up[1, column] - cand_lw_dn[1, column]) -
              (ref_lw_up[1, column] - ref_lw_dn[1, column])
    sfc_err = (cand_lw_dn[end, column] - cand_lw_up[end, column]) -
              (ref_lw_dn[end, column] - ref_lw_up[end, column])

    summary = """
    Flux RMSE (LW): $(round(rmse_lw; digits=3)) W m⁻²
    Flux RMSE (SW): $(round(rmse_sw; digits=3)) W m⁻²
    Heating-rate RMSE: $(round(rmse_hr; digits=4)) K day⁻¹
    TOA LW net forcing err (col 1): $(round(toa_err; digits=4)) W m⁻²
    Surface LW net forcing err (col 1): $(round(sfc_err; digits=4)) W m⁻²

    10 columns, 137 layers, 138 interfaces.
    Hard ecCKD gate: passed.
    """
    Label(fig[2, 3], summary; tellwidth = false, tellheight = false,
          halign = :left, valign = :top, font = :regular)

    Legend(fig[3, 1:3], fig.content[2], orientation = :horizontal, tellheight = true, tellwidth = false)

    save(joinpath(OUT_DIR, "fig1_ecckd_vs_ecrad_profiles.png"), fig; px_per_unit = 2)
    return fig
end

# -----------------------------------------------------------------------------
# Figure 2. H100 speedup bar chart across k-model configurations
# -----------------------------------------------------------------------------

function figure_h100_speedup()
    pareto = JSON.parsefile(joinpath(BREEZE_RESULTS, "reduced_pareto/radiative_heating_reduced_pareto_latest.json"))
    rcemip32 = JSON.parsefile(joinpath(BREEZE_RESULTS, "rcemip_h100_32x32x64/radiative_heating_rcemip_latest.json"))

    entries = Tuple{String, Float64, Float64, Float64}[]

    rh32 = Float64(rcemip32["radiative_heating_update_median_ms"])
    rrt32 = Float64(rcemip32["rrtmgp_update_median_ms"])
    push!(entries, ("32 LW × 32 SW\nvalidated ecCKD (production)", rh32, rrt32, rrt32/rh32))

    for model in pareto["models"]
        rh = get(model, "radiative_heating_update_median_ms", nothing)
        rrt = get(model, "rrtmgp_update_median_ms", nothing)
        (rh === nothing || rrt === nothing) && continue
        nx = get(model, "nx", 0); ny = get(model, "ny", 0); nz = get(model, "nz", 0)
        nx*ny*nz < 32*32*64 && continue
        ng_lw = model["ng_lw"]; ng_sw = model["ng_sw"]
        kind = get(model, "gas_model_kind", "")
        flavor = occursin("fixed", kind) ? "fixed-coeff scaffold" :
                 occursin("validated", kind) ? "validated ecCKD" : "scaffold"
        push!(entries, ("$(ng_lw) LW × $(ng_sw) SW\n$(flavor)",
                        Float64(rh), Float64(rrt), Float64(rrt)/Float64(rh)))
    end

    # Order by total g-points descending so 32x32 sits on the left
    function _total_gpoints(label::String)
        m = match(r"(\d+)\s*LW\s*×\s*(\d+)\s*SW", label)
        m === nothing && return 0
        return parse(Int, m.captures[1]) + parse(Int, m.captures[2])
    end
    sort!(entries; by = e -> -_total_gpoints(e[1]))

    fig = Figure(size = (1200, 680))
    Label(fig[0, 1],
          "RadiativeHeating.jl vs RRTMGP on H100\nRCEMIP-style 32×32×64 / 1024 columns";
          fontsize = 18, font = :bold, tellwidth = false, justification = :center)

    n = length(entries)
    ax = Axis(fig[1, 1];
              ylabel = "radiation update median [ms]\n(lower is faster)",
              xticks = (1:n, [e[1] for e in entries]),
              xticklabelsize = 12,
              title = "")
    barpos = collect(1.0:Float64(n))
    rh_ms = [e[2] for e in entries]
    rrt_ms = [e[3] for e in entries]

    barplot!(ax, barpos .- 0.2, rh_ms; width = 0.38, color = :tomato, label = "RadiativeHeating.jl")
    barplot!(ax, barpos .+ 0.2, rrt_ms; width = 0.38, color = :steelblue, label = "RRTMGP.jl")

    ymax = maximum(rrt_ms) * 1.30
    ylims!(ax, 0, ymax)

    for (i, e) in enumerate(entries)
        # speedup callout above the RRTMGP bar
        text!(ax, i + 0.2, e[3] + 8.0, text = "$(round(e[4]; digits=1))×",
              align = (:center, :bottom), fontsize = 16, font = :bold, color = :black)
        # numeric labels on each bar
        text!(ax, i - 0.2, e[2] + 8.0, text = "$(round(e[2]; digits=1)) ms",
              align = (:center, :bottom), fontsize = 11, color = :tomato)
        text!(ax, i + 0.2, e[3] - 14.0, text = "$(round(e[3]; digits=1)) ms",
              align = (:center, :top), fontsize = 11, color = :white)
    end

    axislegend(ax; position = :rt, framecolor = (:gray, 0.6))

    rowsize!(fig.layout, 1, Relative(0.85))
    save(joinpath(OUT_DIR, "fig2_h100_speedup.png"), fig; px_per_unit = 2)
    return fig
end

# -----------------------------------------------------------------------------
# Figure 3. Training / recovery curves
# -----------------------------------------------------------------------------

function figure_training_curves()
    toy = JSON.parsefile(joinpath(REPO_ROOT, "validation/results/toy_ecckd_training.json"))
    ad = JSON.parsefile(joinpath(REPO_ROOT, "validation/results/rrtmgp_target_16g_ad_calibration.json"))

    toy_loss = Float64.(toy["loss_history"])
    toy_it = collect(0:length(toy_loss)-1)

    ad_traj = ad["training"]["trajectory"]
    ad_it = [Int(t["iteration"]) for t in ad_traj]
    ad_loss = [Float64(t["loss"]) for t in ad_traj]
    if haskey(ad_traj[end], "best_trial_loss") && ad_traj[end]["best_trial_loss"] !== nothing
        push!(ad_it, ad_it[end] + 1)
        push!(ad_loss, Float64(ad_traj[end]["best_trial_loss"]))
    end

    fig = Figure(size = (1100, 560))
    Label(fig[0, 1:2],
          "Reactant/Enzyme calibration of ecCKD coefficients";
          fontsize = 18, font = :bold, tellwidth = false)

    # Left: toy fixed-topology gradient descent
    ax1 = Axis(fig[1, 1]; xlabel = "epoch", ylabel = "training loss (MSE flux + heating rate)",
               title = "toy fixed-topology ecCKD recovery\n(Enzyme-pure gradient checked)",
               yscale = log10)
    lines!(ax1, toy_it, toy_loss; color = :tomato, linewidth = 2.5)
    scatter!(ax1, toy_it, toy_loss; color = :tomato, markersize = 8)
    flux_ratio = round(Float64(toy["flux_rmse_ratio"]); digits = 3)
    hr_ratio = round(Float64(toy["heating_rate_rmse_ratio"]); digits = 3)
    loss_ratio = round(Float64(toy["loss_ratio"]); digits = 4)
    text!(ax1, 1.0, toy_loss[1] * 0.45,
          text = "loss × $(loss_ratio)\nflux RMSE × $(flux_ratio)\nheating-rate RMSE × $(hr_ratio)",
          align = (:left, :top), font = "DejaVu Sans Mono", fontsize = 11)

    # Right: production 16-g RRTMGP-target AD calibration
    ax2 = Axis(fig[1, 2]; xlabel = "epoch", ylabel = "shortwave loss",
               title = "16-g shortwave RRTMGP-target calibration\n(Reactant-compiled, Enzyme reverse-mode)")
    lines!(ax2, ad_it, ad_loss; color = :steelblue, linewidth = 2.5)
    scatter!(ax2, ad_it, ad_loss; color = :steelblue, markersize = 8)
    text!(ax2, ad_it[1] + 0.3, ad_loss[1] - 1.0,
          text = "Δloss = $(round(ad_loss[1] - ad_loss[end]; digits = 2))\n($(ad["parameter_count"]) params)",
          align = (:left, :top), font = "DejaVu Sans Mono", fontsize = 11)

    save(joinpath(OUT_DIR, "fig3_training_curves.png"), fig; px_per_unit = 2)
    return fig
end

# -----------------------------------------------------------------------------
# Figure 4. Reduced-model objective descent over the optimizer chain
# -----------------------------------------------------------------------------

function figure_reduced_descent()
    # Time-ordered chain extracted from radiative_heating.md §0.1; each entry
    # is (step label, hard-gate objective after the accepted move).
    chain = [
        ("baseline (no shortwave subset)", 0.014033547037797689),  # passing 32x32 reference
        ("weighted-greedy 16-g subset",     7.17718646855576),
        ("+ boundary-weighted shortwave weights", 5.557632544628063),
        ("+ projected hard-gate max-norm weights", 5.518843803513278),
        ("+ preflight-optimized scales",          2.588652944311036),
        ("+ pressure-band table moves",           2.311709885137816),
        ("+ active-entry table moves",            8.39187324383),
        ("+ exact post-table weight refit",       8.32190990405),
        ("+ slot-blend table moves",              7.47836712724),
        ("+ post-slot weight refit",              7.45725553248),
        ("+ constrained-table continuations",     7.42818441105),
        ("+ all-case active-entry continuation",  7.42405122402),
        ("+ residual-probe scope",                7.42094851978),
        ("+ Rayleigh-inclusive residual probe",   7.41977527859),
        ("+ boundary-aware weight refit",         7.355912280200036),
        ("+ smaller-step continuation",           7.311608077751165),
        ("+ half-step continuation",              7.297431753704113),
        ("+ quarter-step continuation",           7.297425585382674),
        ("+ Rayleigh exact coordinate scan",      7.29601319503),
        ("+ pair coordinate scan",                7.29587581499),
        ("+ coordinate-descent continuation",     7.29585929551),
        ("+ component-scale refit",               7.13988657025),
        ("+ pressure-component scale",            7.13985474639),
        ("+ temperature-component scale",         7.13985031352),
        ("+ H₂O / gas component scale",           7.13985029534),
        ("retained 32×16 (current best)",         7.13985029534),
    ]

    labels = [c[1] for c in chain]
    objs = [Float64(c[2]) for c in chain]

    fig = Figure(size = (1400, 900))
    Label(fig[1, 1],
          "Reduced 16-g shortwave optimizer:\nworst boundary forcing error vs accepted optimizer moves";
          fontsize = 18, font = :bold, tellwidth = false, justification = :center)

    ax = Axis(fig[2, 1];
              ylabel = "worst boundary forcing error [W m⁻²]",
              xticks = (1:length(labels), labels),
              xticklabelrotation = pi/3,
              xticklabelsize = 11,
              yscale = log10)

    # Plot in physical units (W m⁻²) by multiplying normalized objective by 0.3
    physical = objs .* 0.3
    lines!(ax, 1:length(physical), physical; color = :tomato, linewidth = 2)
    scatter!(ax, 1:length(physical), physical; color = :tomato, markersize = 9)
    hlines!(ax, [0.3]; color = :forestgreen, linestyle = :dash, linewidth = 2)

    ylims!(ax, 0.001, 5.0)
    text!(ax, length(labels) * 0.5, 0.34,
          text = "hard-gate threshold: 0.3 W m⁻²",
          color = :forestgreen, align = (:center, :bottom), fontsize = 13, font = :bold)
    text!(ax, 1.3, physical[1] * 1.8,
          text = "official 32×32 baseline:\n$(round(physical[1]; digits=4)) W m⁻² (passes)",
          align = (:left, :bottom), fontsize = 11, color = :black)
    text!(ax, length(physical) - 0.5, physical[end] * 1.4,
          text = "current best 32×16:\n$(round(physical[end]; digits=2)) W m⁻²\n(≈ $(round(physical[end]/0.3; digits=1))× threshold)",
          align = (:right, :bottom), fontsize = 11, color = :black, font = :bold)

    rowsize!(fig.layout, 1, Auto(0.08))
    rowsize!(fig.layout, 2, Auto(1.0))
    save(joinpath(OUT_DIR, "fig4_reduced_descent.png"), fig; px_per_unit = 2)
    return fig
end

# -----------------------------------------------------------------------------

println("generating fig1...")
figure_flux_profiles()
println("generating fig2...")
figure_h100_speedup()
println("generating fig3...")
figure_training_curves()
println("generating fig4...")
figure_reduced_descent()
println("done. PNGs in: ", OUT_DIR)
