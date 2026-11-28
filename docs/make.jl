# Build the documentation:
#
#     julia --startup-file=no --project=docs docs/make.jl

using Documenter
using DocumenterCitations
using DocumenterVitepress
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
        "Symplectic redistribution" => "obstructions.md",
        "Standard experiments" => "standard.md"
    ],
    "Hamiltonian particle projection" => [
        "Overview" => "projection/index.md",
        "Construction" => "projection/construction.md",
        "Non-separable stages" => "projection/nonseparable.md",
        "Charge-conserving reweighting" => "projection/charge.md",
        "Target weights from a low-rank density" => "projection/lowrank.md",
        "Algorithm" => "projection/algorithm.md",
        "Numerical experiments" => "projection/experiments.md",
        "Structural statements" => "projection/checks.md",
        "Hole filling" => "projection/holes.md",
        "Properties, limitations and higher dimensions" => "projection/limitations.md"
    ],
    "Spline-projection redistribution" => [
        "Overview" => "redistribution/index.md",
        "Positive, moment-compatible spline projection" => "redistribution/projection.md",
        "Symmetric ensemble Hamiltonians" => "redistribution/hamiltonians.md",
        "Symplectic step and parameter solver" => "redistribution/solver.md",
        "Closed loops" => "redistribution/closedloop.md",
        "Algorithm" => "redistribution/algorithm.md",
        "Numerical experiments" => "redistribution/experiments.md",
        "The verification set-up" => "redistribution/setup.md",
        "Hole filling" => "redistribution/holes.md",
        "Properties, limitations and higher dimensions" => "redistribution/limitations.md"
    ],
    "Gauge rezoning and canonical transport" => [
        "Overview" => "remapping/index.md",
        "Velocity rezoning by the gauge step" => "remapping/gauge.md",
        "Canonical transport with consistent weights" => "remapping/transport.md",
        "Algorithm" => "remapping/algorithm.md",
        "Numerical experiments" => "remapping/experiments.md",
        "Clustered markers and a separable density" => "remapping/examples.md",
        "Properties, limitations and higher dimensions" => "remapping/limitations.md"
    ],
    "Sliced redistribution" => [
        "Overview" => "sliced/index.md",
        "The stages and their structure" => "sliced/structure.md",
        "Algorithm" => "sliced/algorithm.md",
        "Numerical experiments" => "sliced/experiments.md",
        "Structural suites" => "sliced/suites.md",
        "Long-time run: variants and repeated remaps" => "sliced/long-time.md",
        "Hole filling" => "sliced/holes.md",
        "Constructions that fail" => "sliced/failures.md",
        "Properties, limitations and higher dimensions" => "sliced/limitations.md"
    ],
    "Comparison and open questions" => "comparison.md",
    "Verification" => [
        "Overview and findings" => "verification.md",
        "Hamiltonian particle projection" => "verification/projection.md",
        "Spline-projection redistribution" => "verification/redistribution.md",
        "Gauge rezoning and canonical transport" => "verification/remapping.md",
        "Sliced redistribution" => "verification/sliced.md"
    ],
    "Library" => [
        "Common layer" => "library.md",
        "Projection" => "projection/library.md",
        "Redistribution" => "redistribution/library.md",
        "Remapping" => "remapping/library.md",
        "SlicedRedistribution" => "sliced/library.md"
    ],
    "Figures and animations" => "figures.md",
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
    # no source or edit links
    remotes = nothing,
    # VitePress builds one site per deploy base: docs/build/1 for a local build, a push to main or
    # a pull request, and docs/build/1, docs/build/2, … for the aliases of a release tag. An
    # animation is written next to its page and embedded as `./<name>.mp4`, which VitePress
    # bundles with the page.
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/michakraus/ParticleRemapping.jl",
        devbranch = "main",
        description = "Symplectic remapping, projection and redistribution of weighted particle-in-cell markers"
    )
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/michakraus/ParticleRemapping.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
)
