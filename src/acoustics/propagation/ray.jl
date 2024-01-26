export Ray
export RayModel

"Equations 3.23-24 of Jensen, et al (2011)."
function eikonal!(du, u, p, s, cel_fcn::Function)
    ∂c_∂x(x, z) = derivative(x -> cel_fcn(x, z), x)
    ∂c_∂z(x, z) = derivative(z -> cel_fcn(x, z), z)

    τ, x, z, ξ, ζ = u
    c = cel_fcn(x, z)
    c² = c^2

    du[1] = dτ_ds = 1 / c
    du[2] = dx_ds = c * ξ
    du[3] = dz_ds = c * ζ
    du[4] = dξ_ds = -∂c_∂x(x, z) / c²
    du[5] = dζ_ds = -∂c_∂z(x, z) / c²
end

struct Ray
    s_max

    x
    z

    function Ray(scen::Scenario, θ₀::Real)
        c = scen.env.ocn.cel.fun

        eikonal_local!(du, u, p, s) = eikonal!(du, u, p, s, c)

        τ₀ = 0.0
        x₀ = scen.x
        z₀ = scen.z
        c₀ = c(x₀, z₀)
        ξ₀ = cos(θ₀) / c₀
        ζ₀ = sin(θ₀) / c₀
        u₀ = [τ₀, x₀, z₀, ξ₀, ζ₀]

        prob = ODEProblem(eikonal_local!, u₀, [0.0, 300e3])
        sol = solve(prob)

        x(s) = sol(s, idxs = 2)
        z(s) = sol(s, idxs = 3)

        s_max = sol.t[end]

        new(s_max, x, z)
    end
end

function MakieCore.convert_arguments(plot_type::MakieCore.PointBased, ray::Ray)
    s = range(0.0, ray.s_max, 301)
    MakieCore.convert_arguments(plot_type, ray.x.(s), ray.z.(s))
end

# function MakieCore.convert_arguments(plot_type::MakeCore.PointBased, rays::AbstractVector{<:Ray})
#     Makie
# end

function default_launch_angles(scen::Scenario)
    atan(1e3 / 10e3) * range(-1, 1, 7) # for Munk profile
end

struct RayModel <: Propagation
    scen::Scenario
    rays::Vector{Ray}

    function RayModel(scen::Scenario, launch_angles::AbstractVector{<:Real} = default_launch_angles(scen))
        rays = [Ray(scen, θ₀) for θ₀ in launch_angles]
        # TODO compute pressure field

        new(scen, rays)
    end
end

function Propagation(::Val{:ray}, scen::Scenario)
    RayModel(scen)
end