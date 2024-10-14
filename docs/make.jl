push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"src/")

using GeMotion
using Documenter

makedocs(
  sitename="GeMotion.jl",
  remotes = nothing,
  pages = Any[
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

# optionally deploy docs to github pages
deploydocs(
    repo = "github.com/lamBOOO/GeMotion.jl.git",
)
