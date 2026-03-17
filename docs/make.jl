using MotifDatabaseSearch
using Documenter

DocMeta.setdocmeta!(MotifDatabaseSearch, :DocTestSetup, :(using MotifDatabaseSearch); recursive=true)

makedocs(;
    modules=[MotifDatabaseSearch],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    sitename="MotifDatabaseSearch.jl",
    format=Documenter.HTML(;
        canonical="https://kchu25.github.io/MotifDatabaseSearch.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/MotifDatabaseSearch.jl",
    devbranch="main",
)
