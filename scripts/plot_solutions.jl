using Serialization
using Plots
using SemiLagrangian

polar_plot_set = load_polars("../polars/Sunodyssey40.pol")
polar_plot = polar_plot_set[4]
polar_plot.boatspeeds .*= 0.05 / maximum(polar_plot.boatspeeds)

σlist     = [0.0, 0.05, 0.1, 0.05, 0.05, 0.05]
driftlist = [0.0, 0.0,  0.0, 0.05, 0.10, 0.15]

# Make directories to hold all of the outputs, if they don't exist
if !isdir("Plots")
    mkdir("Plots")
end
if !isdir("Images")
    mkdir("Images")
end
if !isdir("Movies")
    mkdir("Movies")
end

for (σ, drift) in zip(σlist, driftlist)
    println("Generating for ($σ, $drift)")
    problem_params = ProblemParams(0.1, 1.8, 2.0, σ, drift, 0.0)

    valgrid = deserialize("Grids/RWGS_Gran16Sigma$(σ)Drift$(drift).val")
    dirgrid = deserialize("Grids/RWGS_Gran16Sigma$(σ)Drift$(drift).dir")
    switchgrid = deserialize("Grids/RWGS_Gran16Sigma$(σ)Drift$(drift).swt")

    gr()
    valplot = SemiLagrangian.plot_single_valuegrid(valgrid)
    savefig(valplot, "Plots/Sigma$(σ)Drift$(drift)_value.png")
    switchplot = plot_switchgrid(switchgrid, size=(1000, 1000))
    savefig(switchplot, "Plots/Sigma$(σ)Drift$(drift)_switch.png")

    # Pick a (shared) seed for the trajectories with these parameters
    seed = 1111

    # Traj sims
    all_trajs = SemiLagrangian.PolarTrajectory[]
    for starting_pt in default_starting_pts
        traj = liverun(
            polar_plot, switchgrid, dirgrid, problem_params, starting_pt, 0.01;
            seed=seed
        )
        push!(all_trajs, traj)
    end

    println("\tGenerating plots of trajectories...")
    # If you have PGFPlotsX.jl installed, you can uncomment out this line
    #  to produce prettier plots.
    # pgfplotsx()
    traj_plot = plot_trajectory(all_trajs, problem_params)
    savefig("Images/Trajs_Gran16Sigma$(σ)Drift$(drift).png")

    println("\tAnimating...")
    anim_trajs(all_trajs, problem_params, "Movies/Trajs_Gran16Sigma$(σ)Drift$(drift).mp4"; freq=20)
end