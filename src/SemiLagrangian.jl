module SemiLagrangian

using Parameters
import Random.seed!

include("types.jl")
include("plotting.jl")
include("util.jl")
include("valuesolvers.jl")
include("simulation.jl")

export initialize_valuegrid, compute_switchgrid, single_step, solve_value, solve_value_simplified, solve_value_sweep
export plot_valuegrid, plot_switchgrid, plot_valuegrid_quiver, plot_trajectory
export ProblemParams, deterministic_params, stochastic_params, drift_stochastic_params
export default_grid, finer_grid, extrafine_grid, thefinestgrid
export default_starting_pts, liverun
export paper_speed, ext_paper_speed, rose_speed
export load_polars, plot_polar
export anim_trajs, anim_trajs_makie

end