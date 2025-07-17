using OscillatorOptimization
using Documenter

DocMeta.setdocmeta!(OscillatorOptimization, :DocTestSetup, :(using OscillatorOptimization); recursive=true)

makedocs(;
    modules=[OscillatorOptimization],
    authors="Claude Assistant",
    sitename="OscillatorOptimization.jl",
    format=Documenter.HTML(;
        canonical="https://claude-assistant.github.io/OscillatorOptimization.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/claude-assistant/OscillatorOptimization.jl",
    devbranch="main",
)
