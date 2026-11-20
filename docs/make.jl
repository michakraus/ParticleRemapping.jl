# Build the documentation:
#
#     julia --startup-file=no --project=docs docs/make.jl

using Documenter
using DocumenterCitations
using CairoMakie
using ParticleRemapping

CairoMakie.activate!(; type = "png")
set_theme!(remap_theme())

const BIB = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

const PAGES = [
    "Home" => "index.md",
    "Foundations" => [
        "Background" => "background.md",
        "Canonical structure of weighted markers" => "symplectic.md",
        "What a symplectic redistribution can do" => "obstructions.md"
    ],
    "Gauge rezoning and canonical transport" => [
        "Overview" => "remapping/index.md",
        "Velocity rezoning by the gauge step" => "remapping/gauge.md",
        "Canonical transport with consistent weights" => "remapping/transport.md",
        "Numerical experiments" => "remapping/experiments.md",
        "Combining the two constructions" => "remapping/comparison.md",
        "Properties, limitations and higher dimensions" => "remapping/limitations.md"
    ],
    "Hamiltonian particle projection" => [
        "Overview" => "projection/index.md",
        "Construction" => "projection/construction.md",
        "Non-separable stages" => "projection/nonseparable.md",
        "Charge-conserving reweighting" => "projection/charge.md",
        "Target weights from a low-rank density" => "projection/lowrank.md",
        "Numerical experiments" => "projection/experiments.md",
        "Properties, limitations and higher dimensions" => "projection/limitations.md"
    ],
    "Spline-projection redistribution" => [
        "Overview" => "redistribution/index.md",
        "Positive, moment-compatible spline projection" => "redistribution/projection.md",
        "Symmetric ensemble Hamiltonians" => "redistribution/hamiltonians.md",
        "Symplectic step and parameter solver" => "redistribution/solver.md",
        "Closed loops" => "redistribution/closedloop.md",
        "Numerical experiments" => "redistribution/experiments.md",
        "Properties, limitations and higher dimensions" => "redistribution/limitations.md"
    ],
    "Sliced redistribution" => [
        "Overview" => "sliced/index.md",
        "The stages and their structure" => "sliced/structure.md",
        "Algorithm" => "sliced/algorithm.md",
        "Numerical experiments" => "sliced/experiments.md",
        "What does not work" => "sliced/failures.md",
        "Properties, limitations and higher dimensions" => "sliced/limitations.md"
    ],
    "Comparison and open questions" => "comparison.md",
    "Figures and animations" => "figures.md",
    "Verification" => [
        "Overview and findings" => "verification.md",
        "Gauge rezoning and canonical transport" => "verification/remapping.md",
        "Hamiltonian particle projection" => "verification/projection.md",
        "Spline-projection redistribution" => "verification/redistribution.md",
        "Sliced redistribution" => "verification/sliced.md"
    ],
    "Library" => [
        "Common layer" => "library.md",
        "Remapping" => "remapping/library.md",
        "Projection" => "projection/library.md",
        "Redistribution" => "redistribution/library.md",
        "SlicedRedistribution" => "sliced/library.md"
    ],
    "References" => "references.md"
]

makedocs(;
    modules = [ParticleRemapping],
    authors = "Michael Kraus",
    sitename = "ParticleRemapping.jl",
    pages = PAGES,
    pagesonly = true,
    checkdocs = :exports,
    plugins = [BIB],
    # the package is kept local, without a git remote, so no source or edit links
    remotes = nothing,
    format = Documenter.HTML(;
        # plain URLs: an animation is written next to its page's HTML file and linked by its
        # file name, which a page moved to <page>/index.html would no longer find
        prettyurls = false,
        edit_link = nothing,
        mathengine = Documenter.MathJax3(),
        size_threshold_warn = 400 * 1024,
        size_threshold = 2000 * 1024,
        example_size_threshold = 0,
        assets = String[]
    )
)
