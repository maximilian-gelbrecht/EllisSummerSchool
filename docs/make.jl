using EllisSummerSchool
using Documenter

DocMeta.setdocmeta!(EllisSummerSchool, :DocTestSetup, :(using EllisSummerSchool); recursive=true)

makedocs(;
    modules=[EllisSummerSchool],
    authors="Maximilian Gelbrecht <maximilian.gelbrecht@posteo.de> and contributors",
    sitename="EllisSummerSchool.jl",
    format=Documenter.HTML(;
        canonical="https://maximilian-gelbrecht.github.io/EllisSummerSchool.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/maximilian-gelbrecht/EllisSummerSchool.jl",
    devbranch="main",
)
