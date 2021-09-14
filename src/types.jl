mutable struct ValueGrid
    vgrid   :: Array{Float64, 3}
    Nr      :: Int64
    Nθ      :: Int64
    rbounds :: Tuple{Float64, Float64}
    θbounds :: Tuple{Float64, Float64}
    Δr      :: Float64
    Δθ      :: Float64
end

struct SwitchGrid
    sgrid   :: Array{Int8, 3}
    Nr      :: Int64
    Nθ      :: Int64
    rbounds :: Tuple{Float64, Float64}
    θbounds :: Tuple{Float64, Float64}
    Δr      :: Float64
    Δθ      :: Float64
end

mutable struct DirGrid
    dgrid   :: Array{Float64, 3}
    Nr      :: Int64
    Nθ      :: Int64
    rbounds :: Tuple{Float64, Float64}
    θbounds :: Tuple{Float64, Float64}
    Δr      :: Float64
    Δθ      :: Float64
end

# Parameters defining the problem
struct ProblemParams
    target_radius :: Float64 # radius of the target
    target_dist   :: Float64 # initial distance of target
    switch_cost   :: Float64 # switching cost
    σ             :: Float64 # stochastic diffusion of wind
    drift         :: Float64 # angular drift of wind
    exit_cost     :: Float64 # exit cost
end

# Holds a single polar plot in memory
struct PolarPlot
    angles     :: Vector{Float64}
    boatspeeds :: Vector{Float64}
    windspeed  :: Float64
end

function interpolate_speed(p :: PolarPlot, θ :: Float64) :: Float64
    """ Obtains a boat speed from a polar plot at the given angle
          by linearly interpolating between the available angles.
    """
    θ = mod2pi(θ)
    if θ == 0.0
        return 0.0
    end
    # Linear search for the first angle that is larger than θ
    for (i, pθ) in enumerate(p.angles)
        if θ <= pθ
            frac = (θ - p.angles[i-1]) / (pθ - p.angles[i-1])
            return p[i] * frac + p[i-1] * (1.0 - frac)
        end
    end
    throw(ArgumentError("Argument θ must be in (0, π]"))
end

function resample_polar_plot(p :: PolarPlot, angles :: Vector{Float64}) :: PolarPlot
    """ Given a PolarPlot, returns a new PolarPlot sampled at the angles
         given using linear interpolation between existing angles.
    """
    boatspeeds = [interpolate_speed(p, θ) for θ in angles]
    return PolarPlot(angles, boatspeeds, p.windspeed)
end

function quad_resolution_polar_plot(p :: PolarPlot) :: PolarPlot
    function split(angles)
        new_angles = Vector{Float64}()
        for (i, θ) in enumerate(angles[1:end-1])
            push!(new_angles, θ)
            push!(new_angles, (θ + angles[i+1]) / 2)
        end
        push!(new_angles, angles[end])
        return new_angles
    end

    return resample_polar_plot(p, split(split(p.angles)))
end

function Base.getindex(p :: PolarPlot, i :: Int) :: Float64
    p.boatspeeds[i]
end

function Base.length(p :: PolarPlot)
    length(p.angles)
end

struct PolarTrajectory
    rs :: Vector{Float64}
    θs :: Vector{Float64}
    windθs :: Vector{Float64}
    switch_pts :: Vector{Tuple{Float64, Float64}}
end

function Base.getindex(p :: PolarTrajectory, i :: Int) :: Tuple{Float64, Float64, Float64}
    (p.rs[i], p.θs[i], p.windθs[i])
end

function Base.length(p :: PolarTrajectory)
    length(p.rs)
end

struct CartTrajectory
    xs :: Vector{Float64}
    ys :: Vector{Float64}
    windθs :: Vector{Float64}
    switch_pts :: Vector{Tuple{Float64, Float64}}
end

function Base.convert(CartTrajectory, p :: PolarTrajectory, target_dist)
    xs = @. p.rs * sin(p.θs)
    ys = @. target_dist - p.rs * cos(p.θs)
    switchxs = @. first(p.switch_pts) * sin(last(p.switch_pts))
    switchys = @. target_dist - first(p.switch_pts) * cos(last(p.switch_pts))
    CartTrajectory(xs, ys, p.windθs, collect(zip(switchxs, switchys)))
end

function Base.getindex(p :: CartTrajectory, i :: Int) :: Tuple{Float64, Float64, Float64}
    (p.xs[i], p.ys[i], p.windθs[i])
end

function Base.length(p :: CartTrajectory)
    length(p.xs)
end

# Holds a set of polar plots (from a .pol file) in memory
struct PolarPlotSet
    polars :: Vector{PolarPlot}
    windspeeds :: Vector{Float64}
end

function Base.getindex(p :: PolarPlotSet, i :: Int) :: PolarPlot
    p.polars[i]
end

""" Loads a file containing a set of polar plots (.pol) into a PolarPlotSet.
    Converts angles to radians.
"""
function load_polars(fname :: String) :: PolarPlotSet
    contents = open(readlines, fname)
    num_angles = length(contents) - 1

    windspeed_header = split(contents[1], "\t", keepempty=false)[2:end]
    num_windspeeds = length(windspeed_header)
    windspeeds = map(x -> parse(Float64, x), windspeed_header)

    angles = Array{Float64}(undef, num_angles)
    polardata = Array{Float64}(undef, num_windspeeds, num_angles)

    for i in 1:num_angles
        line_split = split(contents[i+1], "\t", keepempty=false)
        angles[i] = parse(Float64, line_split[1]) * π / 180.0
        polardata[1:end, i] = map(x -> parse(Float64, x), line_split[2:end])
    end

    PolarPlotSet(
        [PolarPlot(angles, polardata[i, 1:end], windspeeds[i]) for i in 1:num_windspeeds],
        windspeeds
    )
end