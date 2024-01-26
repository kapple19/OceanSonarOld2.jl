export Ray
export RayModel

"Equations 3.23-24 of Jensen, et al (2011)."
function eikonal!(du, u, p, s, c::Function)
    ∂c_∂x(x, z) = derivative(x -> c(x, z), x)
    ∂c_∂z(x, z) = derivative(z -> c(x, z), z)

    x, z, ξ, ζ = u

    du[1] = dx_ds = c(x, z) * ξ
    du[2] = dz_ds = c(x, z) * ζ
    du[3] = dξ_ds = -∂c_∂x(x, z) / c(x, z)^2
    du[4] = dζ_ds = -∂c_∂z(x, z) / c(x, z)^2
end

struct Ray
    s_max

    x
    z

    function Ray(scen::Scenario, θ₀::Real)
        c = scen.env.ocn.cel.fun

        eikonal_local!(du, u, p, s) = eikonal!(du, u, p, s, c)

        x₀ = scen.x
        z₀ = scen.z
        c₀ = c(x₀, z₀)
        ξ₀ = cos(θ₀) / c₀
        ζ₀ = sin(θ₀) / c₀
        u₀ = [x₀, z₀, ξ₀, ζ₀]

        prob = ODEProblem(eikonal_local!, u₀, [0.0, 300e3])
        sol = solve(prob)

        x(s) = sol(s, idxs = 1)
        z(s) = sol(s, idxs = 2)

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