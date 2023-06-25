using BiotBrinkmanWithVorticityPaper
using Documenter

DocMeta.setdocmeta!(BiotBrinkmanWithVorticityPaper, :DocTestSetup, :(using BiotBrinkmanWithVorticityPaper); recursive=true)

makedocs(;
    modules=[BiotBrinkmanWithVorticityPaper],
    authors="Ruben Caraballo-Diaz <ruben.caraballo@alumnos.ubiobio.cl>, Chansophea Wathanak In <cwin18@student.monash.edu>, Alberto F. Martin <alberto.f.martin@anu.edu.au>, Ricardo Ruiz-Baier <ricardo.ruizbaier@monash.edu>",
    repo="https://github.com/amartinhuertas/BiotBrinkmanWithVorticityPaper.jl/blob/{commit}{path}#{line}",
    sitename="BiotBrinkmanWithVorticityPaper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://amartinhuertas.github.io/BiotBrinkmanWithVorticityPaper.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amartinhuertas/BiotBrinkmanWithVorticityPaper.jl",
    devbranch="main",
)
