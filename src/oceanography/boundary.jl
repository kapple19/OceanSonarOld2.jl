export Boundary

struct Boundary <: OcnSon
    fun
end

function (bnd::Boundary)(x::Real)
    bty.fun(x)
end