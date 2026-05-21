module ReducedEcCKDMetadataTestHelpers
include(joinpath(@__DIR__, "..", "validation", "reduced_ecckd_accuracy.jl"))
end

@testset "reduced ecCKD model metadata" begin
    helpers = ReducedEcCKDMetadataTestHelpers
    model = helpers.EcCKDTabulatedGasOpticsModel(
        gas_names = (:h2o,),
        pressure_grid = [1000.0, 100000.0],
        temperature_grid = [220.0, 300.0],
        longwave_absorption = fill(0.01, 4, 1, 2, 2),
        shortwave_absorption = fill(0.02, 4, 1, 2, 2),
        shortwave_h2o_absorption = fill(0.03, 4, 2, 2, 1),
        longwave_weights = fill(0.25, 4),
        shortwave_weights = [0.1, 0.2, 0.3, 0.4],
        longwave_source_scale = ones(4),
    )

    reduced = helpers.indexed_tabulated_model(model, [1, 3], [2, 4])
    @test !helpers.uses_full_official_shortwave_weights(reduced)
    @test helpers.uses_full_official_shortwave_weights(model)

    albedo = reshape(collect(1.0:8.0), 4, 2)
    direct_albedo = reshape(collect(21.0:28.0), 4, 2)
    incoming = reshape(collect(11.0:18.0), 4, 2)
    reduced_albedo, reduced_direct_albedo, reduced_incoming =
        helpers.reduced_shortwave_boundary_arrays(albedo, direct_albedo, incoming, reduced)

    @test reduced_albedo == albedo[[2, 4], :]
    @test reduced_direct_albedo == direct_albedo[[2, 4], :]
    @test reduced_incoming == incoming[[2, 4], :]

    groups = helpers.cumulative_weight_groups([0.45, 0.05, 0.25, 0.05, 0.15, 0.05], 3)
    @test length(groups) == 3
    @test all(!isempty, groups)
    @test sort(vcat(groups...)) == collect(1:6)

    moves = (
        (
            component = "static_absorption",
            local_gpoint_index = 1,
            gpoint = 2,
            band = 1,
            pressure_index_start = 1,
            pressure_index_end = 1,
            log_scale = log(2.0),
        ),
        (
            component = "dynamic_h2o",
            local_gpoint_index = 2,
            gpoint = 4,
            band = 2,
            pressure_index_start = 2,
            pressure_index_end = 2,
            log_scale = log(3.0),
        ),
    )
    movable = helpers.indexed_tabulated_model(model, [1, 3], [2, 4])
    helpers.apply_pressure_band_table_moves!(movable, moves)
    @test movable.shortwave_absorption[1, 1, 1, 1] ≈ 0.04
    @test movable.shortwave_absorption[1, 1, 2, 1] ≈ 0.02
    @test movable.shortwave_h2o_absorption[2, 2, 1, 1] ≈ 0.09
    @test movable.shortwave_h2o_absorption[2, 1, 1, 1] ≈ 0.03

    active_moves = (
        (
            component = "static_absorption",
            local_gpoint_index = 1,
            gpoint = 2,
            gas_index = 1,
            pressure_index = 1,
            temperature_index = 2,
            h2o_index = 0,
            log_scale = log(5.0),
        ),
        (
            component = "dynamic_h2o",
            local_gpoint_index = 2,
            gpoint = 4,
            gas_index = 0,
            pressure_index = 2,
            temperature_index = 1,
            h2o_index = 1,
            log_scale = log(7.0),
        ),
    )
    active_movable = helpers.indexed_tabulated_model(model, [1, 3], [2, 4])
    helpers.apply_active_table_entry_moves!(active_movable, active_moves)
    @test active_movable.shortwave_absorption[1, 1, 1, 2] ≈ 0.1
    @test active_movable.shortwave_absorption[1, 1, 1, 1] ≈ 0.02
    @test active_movable.shortwave_h2o_absorption[2, 2, 1, 1] ≈ 0.21
    @test active_movable.shortwave_h2o_absorption[2, 1, 1, 1] ≈ 0.03

    mktemp() do path, io
        write(io, """
        {
          "pressure_band_table_refinement": {
            "accepted_moves": [{
              "component": "static_absorption",
              "local_gpoint_index": 2,
              "gpoint": 4,
              "band": 3,
              "pressure_index_start": 28,
              "pressure_index_end": 40,
              "log_scale": 0.125,
              "scale": 1.1331484530668263
            }],
            "trajectory": []
          }
        }
        """)
        close(io)
        parsed = helpers.latest_preflight_pressure_band_table_moves(; path)
        @test length(parsed) == 1
        @test parsed[1].component == "static_absorption"
        @test parsed[1].local_gpoint_index == 2
        @test parsed[1].pressure_index_start == 28
        @test parsed[1].pressure_index_end == 40
        @test parsed[1].log_scale == 0.125
    end

    mktemp() do path, io
        write(io, """
        {
          "active_table_entry_refinement": {
            "accepted_moves": [{
              "component": "dynamic_h2o",
              "local_gpoint_index": 2,
              "gpoint": 4,
              "gas_index": 0,
              "pressure_index": 28,
              "temperature_index": 3,
              "h2o_index": 2,
              "log_scale": -0.0625,
              "scale": 0.9394130628134758
            }],
            "trajectory": []
          }
        }
        """)
        close(io)
        parsed = helpers.latest_preflight_active_table_entry_moves(; path)
        @test length(parsed) == 1
        @test parsed[1].component == "dynamic_h2o"
        @test parsed[1].local_gpoint_index == 2
        @test parsed[1].pressure_index == 28
        @test parsed[1].temperature_index == 3
        @test parsed[1].h2o_index == 2
        @test parsed[1].log_scale == -0.0625
    end

    mktemp() do path, io
        write(io, """
        {
          "accepted": true,
          "accepted_moves": [{
            "component": "static_absorption",
            "local_gpoint_index": 1,
            "gpoint": 1,
            "gas_index": 1,
            "pressure_index": 2,
            "temperature_index": 3,
            "h2o_index": 0,
            "log_scale": -0.001953125,
            "scale": 0.9980487811074755
          }]
        }
        """)
        close(io)
        parsed = helpers.latest_boundary_base_constrained_table_optimizer_moves(; path)
        @test length(parsed) == 1
        @test parsed[1].component == "static_absorption"
        @test parsed[1].local_gpoint_index == 1
        @test parsed[1].pressure_index == 2
        @test parsed[1].temperature_index == 3
        @test parsed[1].log_scale == -0.001953125
    end

    mktemp() do path, io
        write(io, """
        {
          "accepted": true,
          "accepted_moves": [{
            "component": "static_absorption",
            "local_gpoint_index": 2,
            "gpoint": 4,
            "gas_index": 1,
            "pressure_index": 4,
            "temperature_index": 5,
            "h2o_index": 0,
            "log_scale": 0.00390625,
            "scale": 1.003913894324853
          }]
        }
        """)
        close(io)
        parsed = helpers.latest_boundary_table_continuation_optimizer_moves(; path)
        @test length(parsed) == 1
        @test parsed[1].component == "static_absorption"
        @test parsed[1].local_gpoint_index == 2
        @test parsed[1].gpoint == 4
        @test parsed[1].pressure_index == 4
        @test parsed[1].temperature_index == 5
        @test parsed[1].log_scale == 0.00390625
    end

    mktemp() do path, io
        write(io, """
        {
          "accepted": true,
          "accepted_moves": [{
            "component": "rayleigh",
            "local_gpoint_index": 3,
            "gpoint": 9,
            "gas_index": 0,
            "pressure_index": 0,
            "temperature_index": 0,
            "h2o_index": 0,
            "log_scale": -0.0078125,
            "scale": 0.9922179382602435
          }]
        }
        """)
        close(io)
        parsed = helpers.latest_boundary_table_coordinate_scan_moves(; path)
        @test length(parsed) == 1
        @test parsed[1].component == "rayleigh"
        @test parsed[1].local_gpoint_index == 3
        @test parsed[1].gpoint == 9
        @test parsed[1].log_scale == -0.0078125
    end

    mktemp() do path, io
        write(io, """
        {
          "accepted": true,
          "accepted_moves": [{
            "component": "rayleigh",
            "local_gpoint_index": 6,
            "gpoint": 13,
            "gas_index": 0,
            "pressure_index": 0,
            "temperature_index": 0,
            "h2o_index": 0,
            "log_scale": -0.0009765625,
            "scale": 0.9990239141819757
          }, {
            "component": "rayleigh",
            "local_gpoint_index": 7,
            "gpoint": 14,
            "gas_index": 0,
            "pressure_index": 0,
            "temperature_index": 0,
            "h2o_index": 0,
            "log_scale": -0.0009765625,
            "scale": 0.9990239141819757
          }]
        }
        """)
        close(io)
        parsed = helpers.latest_boundary_table_pair_coordinate_scan_moves(; path)
        @test length(parsed) == 2
        @test all(move -> move.component == "rayleigh", parsed)
        @test parsed[1].local_gpoint_index == 6
        @test parsed[2].local_gpoint_index == 7
        @test parsed[1].log_scale == -0.0009765625
    end

    mktemp() do path, io
        write(io, """
        {
          "accepted": true,
          "accepted_moves": [{
            "component": "static_absorption",
            "local_gpoint_index": 4,
            "gpoint": 10,
            "gas_index": 1,
            "pressure_index": 52,
            "temperature_index": 4,
            "h2o_index": 0,
            "log_scale": -0.0078125,
            "scale": 0.9922179382602435
          }, {
            "component": "rayleigh",
            "local_gpoint_index": 9,
            "gpoint": 21,
            "gas_index": 0,
            "pressure_index": 0,
            "temperature_index": 0,
            "h2o_index": 0,
            "log_scale": -0.00390625,
            "scale": 0.9961013694701175
          }]
        }
        """)
        close(io)
        parsed = helpers.latest_boundary_table_coordinate_descent_moves(; path)
        @test length(parsed) == 2
        @test parsed[1].component == "static_absorption"
        @test parsed[1].pressure_index == 52
        @test parsed[2].component == "rayleigh"
        @test parsed[2].local_gpoint_index == 9
    end

    mktemp() do path, io
        weights = join(string.(fill(1 / 16, 16)), ", ")
        write(io, """
        {
          "accepted": true,
          "initial_objective": 8.197241828881253,
          "final_objective": 8.183370012114846,
          "final_weights": [$weights]
        }
        """)
        close(io)
        parsed = helpers.latest_post_constrained_weight_refit_weights(; path)
        @test parsed !== nothing
        @test parsed.objective == 8.183370012114846
        @test length(parsed.weights) == 16
        @test sum(parsed.weights) ≈ 1.0
    end

    mktemp() do path, io
        weights = join(string.(fill(1 / 16, 16)), ", ")
        write(io, """
        {
          "accepted": false,
          "initial_objective": 8.183370012114846,
          "final_objective": 8.183370012114846,
          "final_weights": [$weights]
        }
        """)
        close(io)
        parsed = helpers.latest_post_constrained_weight_refit_weights(; path)
        @test parsed !== nothing
        @test parsed.objective == 8.183370012114846
    end

    mktemp() do path, io
        weights = join(string.(fill(1 / 16, 16)), ", ")
        write(io, """
        {
          "accepted": false,
          "initial_objective": 8.183370012114846,
          "final_objective": 8.2,
          "final_weights": [$weights]
        }
        """)
        close(io)
        @test helpers.latest_post_constrained_weight_refit_weights(; path) === nothing
    end

    mktemp() do path, io
        weights = join(string.(fill(1 / 16, 16)), ", ")
        write(io, """
        {
          "accepted": true,
          "initial_objective": 7.478367127241938,
          "final_objective": 7.457255532480841,
          "final_weights": [$weights]
        }
        """)
        close(io)
        parsed = helpers.latest_post_slot_weight_refit_weights(; path)
        @test parsed !== nothing
        @test parsed.objective == 7.457255532480841
        @test length(parsed.weights) == 16
        @test sum(parsed.weights) ≈ 1.0
    end

    mktemp() do path, io
        weights = join(string.(fill(1 / 16, 16)), ", ")
        write(io, """
        {
          "accepted": false,
          "initial_objective": 7.478367127241938,
          "final_objective": 7.478367127241938,
          "final_weights": [$weights]
        }
        """)
        close(io)
        @test helpers.latest_post_slot_weight_refit_weights(; path) === nothing
    end
end
