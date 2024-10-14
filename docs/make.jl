push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"src/")

using GenMatFlow
using Documenter

makedocs(
  sitename="GenMatFlow.jl",
  remotes = nothing,
  pages = Any[
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

# optionally deploy docs to github pages
deploydocs(
    repo = "github.com/lamBOOO/GenMatFlow.jl.git",
)
