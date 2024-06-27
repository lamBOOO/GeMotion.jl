using Test

using GFluxx

@testset "Unit tests" begin
  GFluxx.simulate()
end

@testset "Test Examples" begin
  include("../examples/validation-basak/basak.jl")
end
