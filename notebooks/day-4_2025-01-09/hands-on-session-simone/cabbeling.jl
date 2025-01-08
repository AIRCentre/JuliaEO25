# # Cabbeling simulation
#
# Convection driven by density gradients appearing because of the mixing of 
# cold and hot water at the same denisty. This is cause by the nonlinearity of the equation of state of water. 
# We start with a stable fluid, with hot water (7.55ᵒ C) above cold water (1ᵒ C). 
# The fluid is initially at rest because the density is constant in the domain, but mixing of cold and hot water
# generates temperatures that correspond to higher density (the maximum density of fresh water is at 4ᵒ C) and
# the mixed water starts to sink.

using Oceananigans
using Oceananigans.Models: seawater_density
using SeawaterPolynomials: TEOS10EquationOfState

# We start with defining a 2D grid in the x-z plane with 512x256 grid points.
# This is not enough to resolve the Kolmogorov scale so we are in an LES regime

Nx, Ny = 256, 256
grid = RectilinearGrid(CPU(),
                       size = (Nx, Ny),
                       x = (0, 0.5),
                       z = (-0.5, 0.0),
                       topology = (Bounded, Flat, Bounded))

# We use the TEOS10 equation of state to compute the density of seawater.
# This is a nonlinear equation of state that depends on temperature, salinity and pressure.
# We set the salinity to zero and the gravitational acceleration to 9.81 m/s².

equation_of_state = TEOS10EquationOfState(reference_density=1000)
buoyancy = SeawaterBuoyancy(; equation_of_state,
                              constant_salinity=0,
                              gravitational_acceleration=9.81) 

# Since we are in an LES regime we need a LES closure. In this case we use a high-order WENO
# scheme to provide us the implicit diffusion necessary to represent subgrid-scale diffusion in
# the model. Other options would include an explicit closure like
# the Smagorinsky model (`closure = AnisotropicMinimumDissipation()`).

model = NonhydrostaticModel(; grid, 
                              buoyancy, 
                              advection=WENO(order=7),
                              tracers=:T)
            
# Setting the initial conditions to an unstable equilibrium with hot water above cold water.
# and some random perturbations in the velocity field to trigger the instability.

T₁, T₂ = 1, 7.55 # ᵒC
Tᵢ = (x, z) -> z <= -0.25 ? T₁ : T₂
Ξᵢ = (x, z) -> 1e-4 * randn()

using FileIO

download("https://aircentre.github.io/JuliaEO25/img/logo-juliaeo-25.png", "logo-juliaeo-25.png")

img   = FileIO.load("logo-juliaeo-25.png")
alpha = getproperty.(img, :alpha) .|> Float64
alpha = reverse(alpha', dims=2)
alpha = alpha[1:8:end-40, 1:4:end-20]

Tᵢ = [ifelse(alpha[i, j] == 0, T₁, T₂) for i in 1:Nx, j in 1:Ny]

set!(model, T=Tᵢ, u=Ξᵢ, v=Ξᵢ, w=Ξᵢ) 

# We set the time step to 0.2 seconds and the stop time to 100 seconds.

simulation = Simulation(model; Δt = 0.1, stop_time=100)

# We add a callback to print the progress of the simulation every 10 iterations.

function progress(sim) 
    u, v, w = sim.model.velocities
    @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim), ", max(w): ", maximum(abs, w))
end

add_callback!(simulation, progress, IterationInterval(10))

# Let's set an output writer to collect the density and temperature fields every second.

ρ = seawater_density(model)
T = model.tracers.T

output_writer = JLD2OutputWriter(model, (; ρ, T),
                                 filename = "cabbeling",
                                 schedule = TimeInterval(1),
                                 overwrite_existing = true)
                                        
simulation.output_writers[:jld2] = output_writer

# Running the simulation!

run!(simulation)

# ## Visualizing the simulation
#
# Let's use Makie to visualize the density and temperature fields in the x-z plane.

using GLMakie
using Printf

filename = "cabbeling.jld2"

ρt = FieldTimeSeries(filename, "ρ")
Tt = FieldTimeSeries(filename, "T")

Nt = length(ρt)
Nx = size(ρt, 1)

i = Int(Nx / 2)
n = Observable(length(ρt.times)) 
ρ = @lift interior(ρt[$n], i+1:Nx, 1, :)
T = @lift interior(Tt[$n], 1:i,    1, :)
x, y, z = nodes(ρt)

set_theme!(Theme(fontsize=18))
fig = Figure(size=(1000, 800))

ρrange = (minimum(ρt[1]), maximum(ρt))

ax = Axis(fig[1, 2], aspect=1.2, xlabel="x (m)", ylabel="z (m)")
xlims!(ax, 0, 0.5)
ylims!(ax, -0.5, 0)

hm = heatmap!(ax, x[1:i], z, T, colormap=:magma, colorrange=(1.55, 7))
Colorbar(fig[1, 1], hm, label="Temperature (ᵒC)", flipaxis=false)

hm = heatmap!(ax, x[i+1:end], z, ρ, colormap=Makie.Reverse(:grays), colorrange=ρrange)
Colorbar(fig[1, 3], hm, label="Density (kg m⁻³)")

@show title = @sprintf("t = %.3f seconds", ρt.times[n[]])

record(fig, "cabbeling_2d.mp4", 1:Nt, framerate=5) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end