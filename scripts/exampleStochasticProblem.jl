using Plots
# Uncomment out the line below to make prettier (image) plots, but this requires you install
#  PlotlyJS.jl.
# plotlyjs()

# If you instead put `using SemiLagrangian` you won't have to prefix all of the functions
#  with SemiLagrangian, but in this example I'm being explicit about which functions are
#  from the library
import SemiLagrangian

valuegrid = SemiLagrangian.default_grid()

polar_plot_set = SemiLagrangian.load_polars("./polars/Sunodyssey40.pol")
polar_plot = polar_plot_set[4]
# Normalize the boat speeds so that the max speed is 0.05
polar_plot.boatspeeds .*= 0.05 / maximum(polar_plot.boatspeeds)

valuegrid, dirgrid, switchgrid = SemiLagrangian.solve_value_sweep(
    polar_plot, valuegrid, SemiLagrangian.stochastic_params, 0.0, 1e-8; rowwise=true
)

# Create a Vector to hold trajectories as we compute them
# Be sure to feed all trajectories the same seed!
all_trajs = SemiLagrangian.PolarTrajectory[]
for starting_pt in SemiLagrangian.default_starting_pts
    traj = SemiLagrangian.liverun(
        polar_plot, switchgrid, dirgrid, SemiLagrangian.stochastic_params, starting_pt, 0.01;
        seed=12345
    )
    push!(all_trajs, traj)
end

# Plot the trajectories!
traj_plot = SemiLagrangian.plot_trajectory(all_trajs, SemiLagrangian.stochastic_params)
savefig(traj_plot, "stochasticTrajs.png")

# Make an animation!
println("Animating...")
SemiLagrangian.anim_trajs(all_trajs, SemiLagrangian.stochastic_params, "stochasticTrajs.mp4"; freq=20)
