export Celerity

struct Celerity <: OcnSon
	fun::Function
end

function (cel::Celerity)(x::Real, z::Real)
	cel.fun(x, z)
end