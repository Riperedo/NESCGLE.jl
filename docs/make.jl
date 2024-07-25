using Documenter
using NESCGLE

makedocs(
    sitename = "NESCGLE.jl",
    format = Documenter.HTML(),
    modules = [NESCGLE],
    pages = [
        "NESCGLE" => "index.md",
        "Stability Matrix" => "SM.md",
        "Preparation Protocol" => "PP.md",
        "Dynamics" => "Dynamics.md",
        "Liquids theory" => "Equilibrium.md",
        "Utils" => "utils.md"
    ],
    remotes = nothing,
    checkdocs=:exports
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
