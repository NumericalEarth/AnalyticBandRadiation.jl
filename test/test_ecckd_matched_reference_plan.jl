using Test
using JSON

include(joinpath(@__DIR__, "..", "validation", "ecckd_matched_reference_plan.jl"))

@testset "ecCKD matched reference plan" begin
    result = main()
    @test isfile(MATCHED_REFERENCE_PLAN_JSON)
    @test isfile(MATCHED_REFERENCE_PLAN_MD)
    artifact = JSON.parsefile(MATCHED_REFERENCE_PLAN_JSON)
    @test artifact["case"] == "ecckd_matched_reference_plan"
    @test artifact["status"] == "ready_for_published_parity_validation"
    @test length(artifact["required_cases"]) == 16
    @test artifact["missing_case_count"] == 0
    @test all(case["present"] && case["boundary_gpoints_match"]
              for case in artifact["required_cases"])
    @test all(isempty(case["missing_variables"]) for case in artifact["required_cases"])
    @test count(case -> case["all_sky"], artifact["required_cases"]) == 6
    @test Set(case["required_shortwave_gpoints"] for case in artifact["required_cases"]) ==
          Set([32, 64, 96])
    @test Set(case["required_longwave_gpoints"] for case in artifact["required_cases"]) ==
          Set([32, 64])
    for case in artifact["required_cases"]
        overrides = join(case["namelist_overrides"], "\n")
        @test occursin("gas_model_name=\"ECCKD\"", overrides)
        @test occursin("do_save_spectral_flux=true", overrides)
        @test occursin("gas_optics_lw_override_file_name=", overrides)
        @test occursin("gas_optics_sw_override_file_name=", overrides)
        if case["all_sky"]
            @test occursin("sw_solver_name=\"Tripleclouds\"", overrides)
            @test occursin("use_aerosols=true", overrides)
        end
    end
    @test occursin("ecrad_meridian_ecckd_cloudless_noaer_out.nc",
                   join(artifact["existing_ecrad_outputs"], "\n"))
    @test occursin("projection diagnostic", artifact["rationale"])
    md = read(MATCHED_REFERENCE_PLAN_MD, String)
    @test occursin("Status: **ready_for_published_parity_validation**", md)
    @test occursin("published-model accuracy", md)
    @test occursin("Required ecRad Namelist Overrides", md)
end
