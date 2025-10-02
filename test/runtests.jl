using Test

using GeMotion

@testset "Unit tests" begin
  GeMotion.simulate()
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
