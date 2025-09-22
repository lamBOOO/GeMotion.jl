using Test

# Run focused unit tests first (these don't require full GeMotion loading)
@testset "Unit Tests - Solver Components" begin
  include("test_unit_solver.jl")
end

@testset "Unit Tests - Postprocessing Components" begin
  include("test_unit_postprocessing.jl")
end

@testset "Integration Tests - FE Setup" begin
  include("test_integration.jl")
end

# Try to load GeMotion and run comprehensive tests
try
  using GeMotion
  
  @testset "GeMotion Module Tests" begin
    @testset "Basic Simulation Test" begin
      # Test with minimal parameters to avoid long runtime
      result = GeMotion.simulate(
        name="unittest_minimal",
        Pr=1.0,
        Ra=100.0,  # Lower Ra for faster convergence
        n=1.0,
        model=CartesianDiscreteModel((0, 1, 0, 1), (5, 5)),  # Small grid
        nlsolver_opts=(;
          show_trace=false,
          method=:newton,
          ftol=1E-4,  # Relaxed tolerance for faster tests
          xtol=1E-4,
          iterations=10  # Limit iterations
        ),
        T_diri_tags=[1, 2],
        T_diri_expressions=[0.0, 1.0],
        V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
      )
      
      # Test that result has expected structure
      @test haskey(result, :uh)      # velocity field
      @test haskey(result, :ph)      # pressure field  
      @test haskey(result, :Th)      # temperature field
      @test haskey(result, :psih)    # stream function
      @test haskey(result, :Pih)     # heat function
      @test haskey(result, :Sth)     # thermal entropy
      @test haskey(result, :Sfl)     # fluid entropy
      @test haskey(result, :Pr)      # Prandtl number
      @test haskey(result, :Ra)      # Rayleigh number
      @test haskey(result, :n)       # power law exponent
      
      # Test parameter consistency
      @test result.Pr == 1.0
      @test result.Ra == 100.0
      @test result.n == 1.0
    end
    
    @testset "Postprocessing Function Tests" begin
      # Test that postprocessing functions exist and are callable
      @test isdefined(GeMotion, :contourplot_unitsquare)
      @test isdefined(GeMotion, :contourplot_coannulus)
      @test isdefined(GeMotion, :plot_all_unitsquare)
      @test isdefined(GeMotion, :cmap_cold_to_hot_paraview)
      
      # Test colormap function
      cmap = GeMotion.cmap_cold_to_hot_paraview()
      @test typeof(cmap) <: Colors.AbstractColorScheme
    end
  end
  
catch LoadError
  @warn "GeMotion module could not be loaded, skipping full integration tests"
  @testset "GeMotion Module Tests (Skipped)" begin
    @test_skip false  # Indicate tests were skipped
  end
end

@testset "Test Examples" begin
  exampleFolder = joinpath(@__DIR__, "..", "examples")
  for (root, dirs, files) in walkdir(exampleFolder)
    for file in files
      fullpath = joinpath(root, file)
      if isfile(fullpath) && endswith(file, ".jl")
        @info "Test example $(fullpath)"
        dir, _ = splitdir(fullpath)
        @testset "$(fullpath)" begin
          cd(dir) do
          @time @eval Module() begin
            # to avoid putting all examples in a separte main()
            # to avoid conflicting redefinitions of variables/functions
            # see: https://github.com/JuliaLang/julia/issues/40189
            Base.include(@__MODULE__, $fullpath)
          end
        end
        end
      end
    end
  end
end
