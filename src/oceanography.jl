include("oceanography/boundary.jl")
include("oceanography/bathymetry.jl")
include("oceanography/altimetry.jl")

include("oceanography/celerity.jl")
include("oceanography/ocean/celerity.jl")
include("oceanography/seabed/celerity.jl")
include("oceanography/atmosphere/celerity.jl")

include("oceanography/medium.jl")
include("oceanography/ocean.jl")
include("oceanography/seabed.jl")
include("oceanography/atmosphere.jl")

include("oceanography/environment.jl")

export Scenario

struct Scenario <: OcnSon
	env::Environment
	x::Float64
	z::Float64
end

function Scenario(model::Val, x::Real, z::Real)
	env = Environment(model)
	Scenario(env, x, z)
end