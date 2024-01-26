module OceanSonar

using DifferentialEquations: ODEProblem, solve
using ForwardDiff: derivative
import MakieCore

abstract type OcnSon end
abstract type OcnSonOpt end

include("oceanography.jl")
include("acoustics.jl")

end # module OceanSonar