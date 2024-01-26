## Ray Tracing Demo
using OceanSonar
using CairoMakie

scen = Scenario(:munk_profile |> Val, 0.0, 1e3)
prop = Propagation(:ray |> Val, scen)

fig = Figure()
axis = Axis(fig[1, 1])
for ray in prop.rays
    lines!(axis, ray)
end
fig