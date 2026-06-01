using JSON

include(joinpath(@__DIR__, "..", "validation", "ecckd_original_objective_terms.jl"))

@testset "ecCKD original objective term capture" begin
    result = main()
    @test result.status == "objective_terms_captured"
    @test isfile(ORIGINAL_OBJECTIVE_TERMS_JSON)
    @test isfile(ORIGINAL_OBJECTIVE_TERMS_MD)

    artifact = JSON.parsefile(ORIGINAL_OBJECTIVE_TERMS_JSON)
    @test artifact["case"] == "ecckd_original_objective_terms"
    @test artifact["status"] == "objective_terms_captured"
    @test artifact["implementation_status"] == "terms_captured_not_yet_recovered"
    @test artifact["longwave"]["ckd_function"] == "calc_cost_function_ckd_lw"
    @test artifact["shortwave"]["ckd_function"] == "calc_cost_function_ckd_sw"

    longwave_terms = Set(term["name"] for term in artifact["longwave"]["terms"])
    shortwave_terms = Set(term["name"] for term in artifact["shortwave"]["terms"])
    @test "spectral_boundary_flux" in longwave_terms
    @test "spectral_toa_up_flux_20x" in shortwave_terms
    @test "downwelling_only_heating_rate" in shortwave_terms
    @test all(term["present"] for term in artifact["longwave"]["terms"])
    @test all(term["present"] for term in artifact["shortwave"]["terms"])

    lw_sequence = artifact["pass_sequences"]["longwave"]
    sw_sequence = artifact["pass_sequences"]["shortwave"]
    @test lw_sequence["optimize_modes"] == ["relative-base", "relative-ch4", "relative-n2o", "relative-cfc"]
    @test sw_sequence["optimize_modes"] == ["relative-base", "relative-ch4", "relative-n2o"]
    @test occursin("broadband_weight=0.8", lw_sequence["selected_common_options"])
    @test occursin("broadband_weight=0.4", sw_sequence["selected_common_options"])
    @test all(pass["present"] for pass in lw_sequence["passes"])
    @test all(pass["present"] for pass in sw_sequence["passes"])

    md = read(ORIGINAL_OBJECTIVE_TERMS_MD, String)
    @test occursin("ecCKD Original Objective Terms", md)
    @test occursin("spectral_toa_up_flux_20x", md)
    @test occursin("terms_captured_not_yet_recovered", md)
end
