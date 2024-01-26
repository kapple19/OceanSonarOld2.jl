export Ray
export RayModel

function trace!(du, u, p, s, cel_fcn::Function)
    ∂c_∂x(x, z) = derivative(x -> cel_fcn(x, z), x)
    ∂c_∂z(x, z) = derivative(z -> cel_fcn(x, z), z)

    τ, x, z, ξ, ζ, p_re, p_im, q_re, q_im = u
    c = cel_fcn(x, z)
    c² = c^2
    ∂²c_∂x² = derivative(x -> ∂c_∂x(x, z), x)
    ∂²c_∂z² = derivative(z -> ∂c_∂z(x, z), z)
    ∂²c_∂x∂z = derivative(z -> ∂c_∂x(x, z), z)

    # Equation 3.62 of Jensen, et al (2011)
    cnn = c² * (
        ∂²c_∂x² * ζ^2
        - 2∂²c_∂x∂z * ξ * ζ
        + ∂²c_∂z² * ξ^2
    )

    # Equation 3.31 of Jensen, et al (2011)
    du[1] = dτ_ds = 1 / c

    # Equations 3.23-24 of Jensen, et al (2011)
    du[2] = dx_ds = c * ξ
    du[3] = dz_ds = c * ζ
    du[4] = dξ_ds = -∂c_∂x(x, z) / c²
    du[5] = dζ_ds = -∂c_∂z(x, z) / c²

    # Equation 3.58 of Jensen, et al (2011)
    du[6] = dp_re_ds = -q_re * cnn / c^2
    du[7] = dp_im_ds = -q_im * cnn / c^2
    du[8] = dq_re_ds = p_re * c
    du[9] = dq_im_ds = p_im * c
end

struct Ray
    s_max

    x
    z
    A

    function Ray(scen::Scenario, θ₀::Real)
        c = scen.env.ocn.cel.fun

        trace_local!(du, u, p, s) = trace!(du, u, p, s, c)

        # Jensen, et al (2011)
        τ₀ = 0.0 # Eqn 3.32
        x₀ = scen.x # Eqn 3.27
        z₀ = scen.z # Eqn 3.28
        c₀ = c(x₀, z₀)
        ξ₀ = cos(θ₀) / c₀ # Eqn 3.27
        ζ₀ = sin(θ₀) / c₀ # Eqn 3.28
        p₀ = 1 / c₀ # Eqn 3.63
        q₀ = 0.0 # Eqn 3.63
        u₀ = [τ₀, x₀, z₀, ξ₀, ζ₀, real(p₀), imag(p₀), real(q₀), imag(q₀)]

        prob = ODEProblem(trace_local!, u₀, [0.0, 300e3])
        sol = solve(prob)

        x(s) = sol(s, idxs = 2)
        z(s) = sol(s, idxs = 3)
        q(s) = sol(s, idxs = 8) + im * sol(s, idxs = 9)
        c(s) = c(x(s), z(s))

        # Note: The dynamic ray initial conditions (p₀, q₀) determine the ray amplitude equation for A(s).
        # See section 3.3.5.1 of Jensen, et al (2011).
        function A(s)
            √(
                c(s) * cos(θ₀) / (
                    x * c₀ * q(s)
                )
            ) / 4π
        end

        s_max = sol.t[end]

        new(s_max, x, z, A)
    end
end

function MakieCore.convert_arguments(plot_type::MakieCore.PointBased, ray::Ray)
    s = range(0.0, ray.s_max, 301)
    MakieCore.convert_arguments(plot_type, ray.x.(s), ray.z.(s))
end

# function MakieCore.convert_arguments(plot_type::MakeCore.PointBased, rays::AbstractVector{<:Ray})
#     Makie
# end

struct Beam

end

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