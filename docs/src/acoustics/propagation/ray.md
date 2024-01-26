# Acoustic Ray Tracing

```@docs
OceanSonar.eikonal!
```

```@example
using OceanSonar
using CairoMakie

scen = Scenario(:munk_profile |> Val, 0.0, 1e3)
prop = Propagation(:ray |> Val, scen)

fig = Figure()
ax = Axis(fig[1, 1])
for ray in prop.rays
    lines!(ax, ray)
end
fig

save("ray_trace_eg.png", fig) # hide
nothing # hide
```

![Ray trace example](ray_trace_eg.png)
