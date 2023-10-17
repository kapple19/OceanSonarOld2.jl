using Documenter
using OceanSonar

makedocs(
    sitename = "OceanSonar",
    format = Documenter.HTML(),
    modules = [OceanSonar]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
