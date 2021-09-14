#= Generates all of the solution value/direction/switchgrids 
     examined in our paper.
=# 

using Serialization
using Plots
using SemiLagrangian
using Distributed

# Make directories to hold all of the outputs, if they don't exist
if !isdir("Grids")
    mkdir("Grids")
end

σlist     = [0.0, 0.05, 0.1, 0.05, 0.05, 0.05]
driftlist = [0.0, 0.0,  0.0, 0.05, 0.1,  0.15]
granularities = [1, 2, 4, 8, 16, 32]

polar_plot_set = load_polars("./polars/Sunodyssey40.pol")
polar_plot = polar_plot_set[4]
polar_plot.boatspeeds .*= 0.05 / maximum(polar_plot.boatspeeds)

for granularity in granularities
    Δr = SemiLagrangian.default_Δr / granularity
    Δθ = SemiLagrangian.default_Δθ / granularity

    println("Granularity = $granularity")

    # If you launch julia with `julia -p <N>` where N is an integer,
    #  these simulations will be parallelized across N cores.
    @distributed for (i, (σ, drift)) in collect(enumerate(zip(σlist, driftlist)))
        println("σ = $σ, drift = $drift")
        println("\tStarting valuegrid iterations...")

        problem_params = SemiLagrangian.ProblemParams(0.1, 1.8, 2.0, σ, drift, 0.0)

        # Standard Gauss-Seidel run
        println("\tStandard Gauss-Seidel iterations")
        valuegrid = SemiLagrangian.initialize_valuegrid(
            SemiLagrangian.default_rbounds,
            SemiLagrangian.default_θbounds,
            Δr, Δθ,
            SemiLagrangian.deterministic_params  # (σ, drift) is not pulled from here, only grid details
        )
        valuegrid, dirgrid, switchgrid, num_iter = SemiLagrangian.solve_value_sweep(
            polar_plot, valuegrid, problem_params, 0.0, 1e-8; max_iters=2000, rowwise=false
        )

        Serialization.serialize("Grids/GS_Gran$(granularity)Sigma$(σ)Drift$(drift).val", valuegrid)
        Serialization.serialize("Grids/GS_Gran$(granularity)Sigma$(σ)Drift$(drift).dir", dirgrid)
        Serialization.serialize("Grids/GS_Gran$(granularity)Sigma$(σ)Drift$(drift).swt", switchgrid)

        # Rowwise Gauss-Seidel run
        println("\tRowwise Gauss-Seidel iterations")
        valuegrid = SemiLagrangian.initialize_valuegrid(
            SemiLagrangian.default_rbounds,
            SemiLagrangian.default_θbounds,
            Δr, Δθ,
            SemiLagrangian.deterministic_params  # (σ, drift) is not pulled from here, only grid details
        )
        valuegrid, dirgrid, switchgrid, num_iter = SemiLagrangian.solve_value_sweep(
            polar_plot, valuegrid, problem_params, 0.0, 1e-8; max_iters=2000, rowwise=true
        )

        Serialization.serialize("Grids/RWGS_Gran$(granularity)Sigma$(σ)Drift$(drift).val", valuegrid)
        Serialization.serialize("Grids/RWGS_Gran$(granularity)Sigma$(σ)Drift$(drift).dir", dirgrid)
        Serialization.serialize("Grids/RWGS_Gran$(granularity)Sigma$(σ)Drift$(drift).swt", switchgrid)
    end
end
