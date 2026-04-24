using Test
using AnalyticBandRadiation
using Dates

@testset "AnalyticBandRadiation" begin
    include("test_planck.jl")
    include("test_absorption.jl")
    include("test_williams_longwave.jl")
    include("test_shortwave.jl")
    include("test_zenith.jl")
end
