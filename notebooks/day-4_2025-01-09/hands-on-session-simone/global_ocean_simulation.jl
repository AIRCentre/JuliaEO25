# # Realistic Ocean Simulations in pure Julia
#
# In this tutorial we will build a global ocean simulation using the ClimaOcean and Oceananigans Julia packages.
# The simulation will be at a nominal four-degree with parameterizations for mesoscale eddies and vertical mixing.
# We will set up the grid, download and interpolate the bathymetry, 
# download forcing files, and finally set up and run the simulation

# ### ***Bonus***
# At the end of the tutorial we will change the resolution to allow the spontaneous generation of eddies and 
# remove the eddy parameterization to see some beautiful ocean turbulence develop! (make sure to have GPU access
# for this step!)

# ## Required software
#
# The tutorial is quite computationally expensive, therefore, if you have 
# access, it is recommended to run the tutorial on one GPU. 
# However, for the purpose of understanding how the library works, 
# a 4 degree global ocean model runs comfortably on a laptop.

# ## Packages:

# Let's start by importing the necessary packages, these are:
# - ClimaOcean: the ocean model
# - Oceananigans: the fluid dynamics library doing the heavy lifting
# - OrthogonalSphericalShellGrids: contains the constructor for the Tripolar grid we will use in the tutorial
# - Printf: always useful for spitting out output
# - CairoMakie: visualization package to visualize the results

using Pkg
pkg"add Oceananigans#main, ClimaOcean#main, OrthogonalSphericalShellGrids, CairoMakie"
using ClimaOcean
using Oceananigans
using Oceananigans.Units
using OrthogonalSphericalShellGrids
using CFTime
using Dates
using Printf
using CairoMakie

# # Building a Global Ocean domain

# We will start building a global ocean in steps: 
# - (1) specifying an **Architecture**
# - (2) choosing a **Vertical coordinate** 
# - (2) Building a **Grid** 
# - (3) Downloading an interpolating a **Bathymetry**

# ### Architectures

# using an architecture is easy...
# it is possible to choose between:
# - (1) CPU 
# - (2) GPU 
# - (3) Distributed 

# In this case, we select `CPU()`, which is always the right choice to start prototyping. Building the script on CPUs is better to catch potential bugs in the script. 
# Once we know that everything works correctly, we can just change the following line to `arch = GPU()`, and voilá, the simulation runs on GPUs

arch = CPU()

# ### Vertical coordinates

# Oceananigans currently supports only $z$ and $z^\star$ coordinates.
# ClimaOcean provides a simple utility to build a simple exponential vertical coordinate

depth = 5000meters
Nz    = 10
h     = 3 # e-folding length of the exponential

r_faces = ClimaOcean.exponential_z_faces(; Nz, h, depth)

# To use z-star coordinates we can do `z_faces = ZStarVerticalCoordinate(r_faces)`, 
# but we will stick with z coordinates for now.

z_faces = r_faces 

# ## Building a grid

# ClimaOcean allows building ocean simulations on three different grid types:
# - `RectilinearGrid`s which represent a _box_ or a Cartesian domain
# - `LatitudeLongitudeGrid`s, which discretizes the sphere along latitude and longitude lines
# - `OrthogonalSphericalShellGrid`s which discretize the sphere with two-dimensional coordinates that do not need to follow latitude and longitude lines. 
#    The only constraint is that the grid must be locally orthogonal in the horizontal direction.
#
# `LatitudeLongitudeGrid`s are the easiest grids to work with since coordinates are one-dimensional and globally orthogonal.
# (i.e. latitude depends only on the `j`-index and longitude depends only on the `i`-index)
# However, `LatitudeLongitudeGrid`s have the problem of the zonal spacing approaching zero as we move to the poles. 
#
# For this reason, to represent a global ocean we use a specific type of `OrthogonalSphericalShellGrid`, 
# called `TripolarGrid` that discretizes the sphere as a set of perpendicular ellipses and hyperbolae.
#
# Let's build a coarse `TripolarGrid` (about 4 degree resolution). 
# We pass to the grid, the architecture, the floating point precision, the size of the grid, and the vertical coordinate.

Nx = 256 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 128 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz   = length(r_faces) - 1
grid = TripolarGrid(arch, Float64; size=(Nx, Ny, Nz), z=z_faces)

# #### Let's visualize the grid using CairoMakie

cartesian_nodes, _ = OrthogonalSphericalShellGrids.get_cartesian_nodes_and_vertices(grid, Face(), Face(), Center())
xF, yF, zF = cartesian_nodes

cartesian_nodes, _ = OrthogonalSphericalShellGrids.get_cartesian_nodes_and_vertices(grid, Center(), Center(), Center())
xC, yC, zC = cartesian_nodes

φ = Oceananigans.Grids.φnodes(grid, Center(), Center())
λ = Oceananigans.Grids.λnodes(grid, Center(), Center())

fig = Figure(size=(1200, 300))
axS = Axis3(fig[1, 1]; aspect=(1, 1, 1), azimuth = 7, height=300)
axC = Axis(fig[1, 2])
surface!(axS, xC .* 0.99, yC .* 0.99, zC .* 0.99)
wireframe!(axS, xC, yC, zC, color = :blue)    

[vlines!(axC, i, color = :green, linewidth=0.5) for i in 1:5:Nx]
[hlines!(axC, j, color = :green, linewidth=0.5) for j in 1:5:Ny]
contour!(axC, φ, levels=-80:10:90, color=:blue, linewidth=0.8)
contour!(axC, λ, levels=0:10:360,  color=:blue, linewidth=0.8)
xlims!(axC, 1, Nx-1)
ylims!(axC, 1, Ny-1)

display(fig)

# ## Adding a bathymetry to the grid
#
# ClimaOcean provides a nifty utility to regrid the bathymetry over the grid, the `regrid_bathymetry` function.
# By default ClimaOcean downloads the ETOPO22 bathymetry at 1/60ᵒ resolution (459 MB) from the NOAA servers.
# However, since the servers are quite busy, I have uploaded a lower resolution version file to dropbox.
# !!! NOTE: This will download the ETOPO22 bathymetry, so make sure that you have an internet connection

url = "https://www.dropbox.com/scl/fi/zy1cu64ybn93l67rjgiq0/Downsampled_ETOPO_2022.nc?rlkey=5upqvoxrnljj205amqf663vcw&st=ou8b32tt&dl=0"
filename = isfile("Downsampled_ETOPO_2022.nc") ? "Downsampled_ETOPO_2022.nc" : download(url, "Downsampled_ETOPO_2022.nc")
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, filename, dir="./")

fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1])
heatmap!(ax, interior(bottom_height, :, :, 1), colormap=:deep)
display(fig)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

# # Configuring an Ocean model
#
# ### Numerical details
#
# Oceananigans allows several numerical options. 
# We use a WENO schemes for the advection of momentum and 
# a centered scheme for tracer advection, to avoid implicit diapycnal diffusion of tracers.
# Stability in the momentum field is ensured by the WENO method. For the tracer field, since the centered
# scheme is dispersive, we need to add some explicit diffusion to avoid numerical instabilities.

momentum_advection = WENOVectorInvariant(order=3) 
tracer_advection   = Centered()

free_surface = SplitExplicitFreeSurface(grid; substeps=30) 

# ### Physical parameterizations
#
# We add an advective GM parameterization for mesoscale eddies and a very simple vertical mixing scheme.
# All the closures require passing also the desired floating point precision of the model

using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, 
                                       ExplicitTimeDiscretization, 
                                       DiffusiveFormulation

eddy_closure    = HorizontalScalarDiffusivity(ν=5e3, κ=1e3)
vertical_mixing = ConvectiveAdjustmentVerticalDiffusivity(Float64;
                                                          convective_κz=5, 
                                                          convective_νz=5, 
                                                          background_κz=3e-5, 
                                                          background_νz=1e-3)

closure = (eddy_closure, vertical_mixing) 

# ### Adding a restoring term
#
# Since we do not have a sea ice model, we strongly restore temperature and salinity 
# at the surface to ECCO climatology in the Polar regions 
# To do this we need to be able to download ECCO data. This is done automatically by `ClimaOcean` 
# provided that you have credentials to access the ECCO data. To do so, follow the instructions detailed in
# **https://github.com/CliMA/ClimaOcean.jl/blob/main/src/DataWrangling/ECCO/README.md**.

start = DateTimeProlepticGregorian(1993, 1, 1)
stop  = DateTimeProlepticGregorian(1993, 3, 1)
dates = range(start, stop; step=Month(1))

temperature = ECCOMetadata(:temperature; dates, dir="./")
salinity    = ECCOMetadata(:salinity;    dates, dir="./")

mask  = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(-100, 0))

FT = nothing # ECCORestoring(temperature, grid; mask, rate=1/2days)
FS = nothing # ECCORestoring(salinity,    grid; mask, rate=1/2days)

forcing = NamedTuple() # (T=FT, S=FS)

# ### Building the ocean simulation
# 
# ClimaOcean provides a utility to build an ocean simulation with all the necessary components. 
# The function `ocean_simulation` returns a `Simulation` object of a `HydrostaticFreeSurfaceModel` that has 
# all the necessary components (BC, drag, etc) to run a global ocean simulation.

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         closure, 
                         forcing, 
                         free_surface)

# ### Initialize our Ocean
#
# We use ECCO climatology to initialize the temperature and salinity fields. 
# We can use the metadata we defined earlier to set the initial conditions. 

#=
u_velocity   = ECCOMetadata(:u_velocity,   dates=dates[1], dir="./")
v_velocity   = ECCOMetadata(:v_velocity,   dates=dates[1], dir="./")
surf_height  = ECCOMetadata(:free_surface, dates=dates[1], dir="./")
=#
set!(ocean.model, T=temperature[1], S=salinity[1]) 

# ### Visualizing the initial conditions
# 
# Let's visualize the surface of our initialized ocean model

fig = Figure(size = (2000, 800))
axT = Axis(fig[1, 1])
axS = Axis(fig[1, 2])

hmT = heatmap!(axT, interior(ocean.model.tracers.T, :, :, grid.Nz), colorrange=(-1, 30), colormap=:magma)
hmS = heatmap!(axS, interior(ocean.model.tracers.S, :, :, grid.Nz), colorrange=(28, 38), colormap=:haline)
Colorbar(fig[0, 1], hmT, vertical=false, label="Surface temperature ᵒC")
Colorbar(fig[0, 2], hmS, vertical=false, label="Surface salinity psu")

display(fig)

# # Adding an atmosphere
#
# ClimaOcean is a prototype for a coupled earth system model. 
# It couples an atmosphere to an ocean and computes air-sea fluxes using bulk formulae.
# At the moment, ClimaOcean provides a utility to download the JRA55 atmospheric reanalysis
# and use it as a prescribed atmosphere.
#
# !!! NOTE: This will download the JRA55 atmospheric reanalysis, so make sure that you have an internet connection (and enough disk space)
#
# We use an idealized atmosphere for this tutorial to avoid downloading the JRA55 data (~15GB).

atmos_grid  = LatitudeLongitudeGrid(arch, Float32; size=(320, 200), latitude=(-90, 90), longitude=(0, 360), topology=(Periodic, Bounded, Flat))
atmos_times = range(0, 360Oceananigans.Units.days, length=10)
atmosphere  = PrescribedAtmosphere(atmos_grid, atmos_times)

# We set up the atmosphere with an idealized temperature and wind speed that do not change in time.

Tₐ(λ, φ) = 255 + cosd(φ) * 80
uₐ(λ, φ) = 10 * sind(2φ)^2
qₐ(λ, φ) = cosd(φ) * 3e-2

for t in eachindex(atmos_times)
    set!(atmosphere.tracers.T[t],    Tₐ)
    set!(atmosphere.tracers.q[t],    qₐ)
    set!(atmosphere.velocities.u[t], uₐ)

    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.tracers.T[t])
    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.tracers.q[t])
    Oceananigans.BoundaryConditions.fill_halo_regions!(atmosphere.velocities.u[t])
end

# If we had a realistic atmosphere we would add radiative properties, however, since we do not have 
# downwelling radiation, we set it to nothing

radiation = nothing # Radiation(ocean_albedo = LatitudeDependentAlbedo())

# ### Coupling the atmosphere to the ocean
#
# The `OceanSeaIceModel` is an `AbstractModel` defined in ClimaOcean that couples an ocean to an atmosphere and a sea ice component.
# For the moment, the sea-ice component is not implemented, so we will only couple the ocean to the atmosphere.
# Instead of the sea ice model, we limit the temperature of the ocean to the freezing temperature.

sea_ice = ClimaOcean.FreezingLimitedOceanTemperature()
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

# ### Building the simulation
#
# We build the simulation out of the `earth_model` as we would do for any other Oceananigans model.
# We start with a smallish time-step (5 minutes) and run only for 2 days to dissipate initialization shocks.

earth = Simulation(earth_model; Δt=5minutes, stop_time=2days)

# ### Adding some diagnostics
#
# We add a callback to save surface fields as well as surface fluxes, every 6 hours

u, v, _ = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
s = sqrt(u^2 + v^2)

η = ocean.model.free_surface.η 

earth.output_writers[:surface_tracers] = JLD2OutputWriter(ocean.model, (; T, S, s),
                                                          schedule = TimeInterval(12hours),
                                                          indices = (:, :, grid.Nz),
                                                          overwrite_existing = true,
                                                          filename = "surface_fields.jld2")


earth.output_writers[:free_surface] = JLD2OutputWriter(ocean.model, (; η),
                                                       schedule = TimeInterval(12hours),
                                                       overwrite_existing = true,
                                                       filename = "free_surface.jld2")

Q  = earth.model.fluxes.total.ocean.heat
τx = earth.model.fluxes.total.ocean.momentum.u
τy = earth.model.fluxes.total.ocean.momentum.v
PE = earth.model.fluxes.total.ocean.tracers.S

earth.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, (; Q, τx, τy, PE),
                                                 schedule = TimeInterval(12hours),
                                                 overwrite_existing = true,
                                                 filename = "surface_fluxes.jld2")

# Also, we add a callback to print a message about how the simulation is going

wall_time = [time_ns()]

function progress(earth)
    clock   = earth.model.clock

    maxu = maximum(abs, u)
    maxv = maximum(abs, v)
    maxT = maximum(T)
    minS = minimum(S)
    
    @info @sprintf("Iteration: %d, time: %s, wall_time: %s, max(|u|, |v|): %.2e %.2e max(T): %.2e, min(S): %.2e\n",
                   clock.iteration, prettytime(clock.time), prettytime(1e-9 * (time_ns() - wall_time[1])), maxu, maxv, maxT, minS)

    wall_time[1] = time_ns()
end

add_callback!(earth, progress, IterationInterval(10))

# ### Running the simulation
#
# quite simply

run!(earth)

earth.Δt = 30minutes
earth.stop_time = 100days

run!(earth)

# ## Visualizing the results
#
# We can visualize the results using CairoMakie. We record a video of surface variables and fluxes.
# To load the data we can use Oceananigans' `FieldTimeSeries` object.

using JLD2
using Oceananigans.Grids: halo_size

file  = jldopen("free_surface.jld2")
iters = keys(file["timeseries/t"]) 

Hx, Hy, _ = halo_size(η.grid)
T  = FieldTimeSeries("surface_fields.jld2", "T")
S  = FieldTimeSeries("surface_fields.jld2", "S")
s  = FieldTimeSeries("surface_fields.jld2", "s")

n  = Observable(1)
Tn = @lift(interior(T[$n], :, :, 1))
Sn = @lift(interior(S[$n], :, :, 1))
sn = @lift(interior(s[$n], :, :, 1))
ηn = @lift(file["timeseries/η/" * iters[$n]][Hx+1:end-Hx, Hy+1:end-Hy, 1])

fig = Figure(size = (1800, 800))
axT = Axis(fig[1, 1], title="Surface temperature ᵒC")
axS = Axis(fig[1, 2], title="Surface salinity psu")
axs = Axis(fig[2, 1], title="Surface speed ms⁻¹")
axη = Axis(fig[2, 2], title="Sea surface height m")

λ, φ, z = nodes(T[1])

hmT = heatmap!(axT, Tn, colormap=:magma,  colorrange=(-1, 30))
hmS = heatmap!(axS, Sn, colormap=:haline, colorrange=(25, 40))
hms = heatmap!(axs, sn, colormap=:deep,   colorrange=( 0, 0.8))
hmη = heatmap!(axη, ηn, colormap=:bwr,    colorrange=(-1, 1))

CairoMakie.record(fig, "surface_fields.mp4", 1:length(T.times); framerate=5) do i 
    @info "doing $i of $(length(T.times))"
    n[] = i
end

# let's also visualize the surface fluxes that force the model

Q  = FieldTimeSeries("surface_fluxes.jld2", "Q")
τx = FieldTimeSeries("surface_fluxes.jld2", "τx")
τy = FieldTimeSeries("surface_fluxes.jld2", "τy")
PE = FieldTimeSeries("surface_fluxes.jld2", "PE")

Qn  = @lift(interior(Q[$n],  :, :, 1))
τxn = @lift(interior(τx[$n], :, :, 1))
τyn = @lift(interior(τy[$n], :, :, 1))
PEn = @lift(interior(PE[$n], :, :, 1))

fig  = Figure(size = (1800, 800))
axQ  = Axis(fig[1, 1], title="Net heat flux Wm⁻²")
axPE = Axis(fig[1, 2], title="Net salt flux psu m s⁻¹")
axτx = Axis(fig[2, 1], title="Zonal wind stress Nm⁻²")
axτy = Axis(fig[2, 2], title="Meridional wind stress Nm⁻²")

hmQ  = heatmap!(axQ,  Qn,  colormap=:magma,   colorrange=(-800,  800))
hmPE = heatmap!(axPE, PEn, colormap=:haline,  colorrange=(-1e-5, 5e-5))
hmτx = heatmap!(axτx, τxn, colormap=:balance, colorrange=(-5e-4, 5e-4))
hmτy = heatmap!(axτy, τyn, colormap=:balance, colorrange=(-5e-4, 5e-4))

CairoMakie.record(fig, "surface_fluxes.mp4", 1:length(Q.times); framerate=5) do i 
    @info "doing $i of $(length(Q.times))"
    n[] = i
end

# Let's visualize the internal structure of temperature and salinity
x, y, z = nodes(ocean.model.tracers.T)

fig = Figure(size = (1200, 400))
axT = Axis(fig[1, 1], title="Internal temperature structure ᵒC")
axS = Axis(fig[1, 2], title="Internal salinity structure psu")

contourf!(axT, 1:48, z, interior(mean(ocean.model.tracers.T, dims=1), 1, :, :), colormap=:magma)
contourf!(axS, 1:48, z, interior(mean(ocean.model.tracers.S, dims=1), 1, :, :), colormap=:haline)

# # Running a high-resolution simulation
#
# What are the steps to modify the above script to run an eddying (quarter degree) simulation?