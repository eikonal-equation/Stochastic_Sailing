using Plots, Images, Parameters
import GLMakie

function plot_valuegrid(valuegrid::ValueGrid; proj=:polar)
    """ Generates an image visualizing the value function grid.
    """
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    
    θticks = ["0", "\\pi/2", "\\pi", "3\\pi/2", "2\\pi"]
    ypixels = ceil.(Int, [0, π/2, π, 3π/2, 2π] ./ Δθ)
    xpixels = 0:floor(Int, Nr/4):Nr
    rticks = xpixels .* Δr
    xlabel = "r"
    ylabel = "\\theta"

    rs = LinRange(rbounds..., Nr)
    θs = LinRange(0, 2π, Nθ)

    valmap1 = vgrid[:, :, 1]
    valmap2 = vgrid[:, :, 2]
    if proj == :polar
        valmap1 = transpose(valmap1)
        valmap2 = transpose(valmap2)
        xlabel = ""
        ylabel = ""
        xticks = ([], [])
        yticks = ([Nr], [string(rbounds[2])])
    end

    grid1 = heatmap(θs, rs, valmap1,
                 xticks=xticks, yticks=yticks, title="q=1",
                 xlabel=xlabel, ylabel=ylabel,
                 color=:inferno, proj=proj)
    grid2 = heatmap(θs, rs, valmap2,
                 xticks=xticks, yticks=yticks, title="q=2",
                 xlabel=xlabel, ylabel=ylabel,
                 color=:inferno, proj=proj)
    plot(grid1, grid2, layout=@layout [a b])
end

function plot_single_valuegrid(valuegrid::ValueGrid; proj=:polar)
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    
    θticks = ["0", "\\pi/2", "\\pi", "3\\pi/2", "2\\pi"]
    ypixels = ceil.(Int, [0, π/2, π, 3π/2, 2π] ./ Δθ)
    xpixels = 0:floor(Int, Nr/4):Nr
    rticks = xpixels .* Δr
    xlabel = "r"
    ylabel = "\\theta"

    valmap1 = vgrid[:, :, 1]
    valmap2 = vgrid[:, :, 2]
    if proj == :polar
        valmap1 = transpose(valmap1)
        valmap2 = transpose(valmap2)
        xlabel = ""
        ylabel = ""
        xticks = ([], [])
        yticks = ([Nr], [string(rbounds[2])])
    end

    heatmap(valmap1,
             xticks=xticks, yticks=yticks, title="q=1",
             xlabel=xlabel, ylabel=ylabel,
             color=:inferno, proj=proj, size=(600, 600))
end

function plot_switchgrid(switchgrid::SwitchGrid; proj=:polar, target_radius=0.1, kwargs...)
    """ Generates an image visualizing the switching sets.
    """
    @unpack sgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = switchgrid
    xlabel = "r"
    ylabel = "\\theta"

    # 1 where q = 1 switches, -1 where q = 2 switches
    switchmap = sgrid[:, :, 1] - sgrid[:, :, 2]

    xpixels = 0:floor(Int, Nr/6):Nr
    rticks = xpixels .* Δr
    θticks = ["0", "\\pi/2", "\\pi", "3\\pi/2", "2\\pi"]
    ypixels = ceil.(Int, [0, π/2, π, 3π/2, 2π] ./ Δθ)
    xticks = (xpixels, rticks)
    yticks = (ypixels, θticks)

    rs = LinRange(rbounds..., Nr)
    θs = LinRange(0, 2π, Nθ)

    if proj == :polar
        switchmap = transpose(switchmap)
        xlabel = ""
        ylabel = ""
        xticks = ([], [])
        yticks = ([Nr], [string(rbounds[2])])
    end
    
    heatmap(θs, rs, switchmap, xlabel=xlabel, ylabel=ylabel, color=:coolwarm,
            colorbar=:none, xticks=xticks, yticks=yticks; proj=proj, kwargs...)
    if proj == :polar
        num_pts = 100
        plot!(range(0, 2π, length=num_pts), fill(target_radius, num_pts), lw=2, color=:black, label="")
    end
    plot!(proj=proj, xtickfontsize=20)
end

function plot_polar(polar::PolarPlot; kwargs...)
    # Mirror data from [0, π] to [π, 2π]
    angles = [polar.angles; reverse(2π .- polar.angles)]
    speeds = [polar.boatspeeds; reverse(polar.boatspeeds)]
    plot(angles, speeds, proj=:polar, m=2, label=""; kwargs...)
end

function circle(target_dist, radius, num_pts)
    angles = range(0.0, stop=2π, length=num_pts)
    return radius .* [cos.(angles) sin.(angles)] .+ [0 target_dist]
end

function plot_trajectory(traj::CartTrajectory, params::ProblemParams; num_quivers=0, lw=4, m=5, kwargs...)
    @unpack xs, ys, windθs, switch_pts = traj
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params

    circle_pts = circle(target_dist, target_radius, 100)
    plot(circle_pts[:, 1], circle_pts[:, 2],
         xlims = (-target_dist, target_dist),
         ylims = (-0.2, target_dist + target_radius),
         label="",
         color=:black,
         lw=3)

    plot!(xs, ys, color=:red, label="", lw=lw)
    if num_quivers != 0
        pt_step = floor(Int, length(xs) / (num_quivers + 1))
        # For a subset of the points, plot wind arrows
        windx = 0.3 .* sin.(windθs)[pt_step:pt_step:end-1]
        windy = -0.3 .* cos.(windθs)[pt_step:pt_step:end-1]
        quiver!(xs[pt_step:pt_step:end-1], ys[pt_step:pt_step:end-1], quiver=(windx, windy), kwargs...)
    end

    # Plot switch points
    scatter!(first.(traj.switch_pts), last.(traj.switch_pts), color=:blue, m=m)
end

function plot_trajectory(polar_traj::PolarTrajectory, params::ProblemParams; num_quivers=0, lw=4, m=5, kwargs...)
    traj = convert(CartTrajectory, polar_traj, params.target_dist)
    plot_trajectory(traj, params; num_quivers=num_quivers, lw=lw, m=m, kwargs...)
end

function plot_trajectory(trajs::Vector{CartTrajectory}, params::ProblemParams; num_quivers=0, lw=4, m=5, kwargs...)
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params
    circle_pts = circle(target_dist, target_radius, 100)
    plot(circle_pts[:, 1], circle_pts[:, 2],
         xlims = (-target_dist, target_dist),
         ylims = (-0.2, target_dist + target_radius + 0.05),
         label="",
         color=:black,
         lw=3,
         xtickfontsize=14,
         ytickfontsize=14)
    for traj in trajs
        @unpack xs, ys, windθs, switch_pts = traj
        plot!(xs, ys, color=:red, label="", lw=lw)

        if num_quivers != 0
            pt_step = floor(Int, length(xs) / (num_quivers + 1))
            # For a subset of the points, plot wind arrows
            windx = 0.3 .* sin.(windθs)[pt_step:pt_step:end-1]
            windy = -0.3 .* cos.(windθs)[pt_step:pt_step:end-1]
            quiver!(xs[pt_step:pt_step:end-1], ys[pt_step:pt_step:end-1], quiver=(windx, windy))
        end
    end

    # Plot all switch points afterwards
    for traj in trajs
        scatter!(first.(traj.switch_pts), last.(traj.switch_pts), color=:blue, m=m, label="", kwargs...)
    end
    plot!()
end

function plot_trajectory(polar_trajs::Vector{PolarTrajectory}, params::ProblemParams; num_quivers=0, lw=4, m=5, kwargs...)
    trajs = convert.(CartTrajectory, polar_trajs, params.target_dist)
    plot_trajectory(trajs, params; num_quivers=num_quivers, lw=lw, m=m, kwargs...)
end

function anim_trajs(trajs::Vector{CartTrajectory}, params::ProblemParams, fname;
                    freq=1, markersize=12, linewidth=4, arrowsize=15, height=600, width=600)
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params
    circle_pts = circle(target_dist, target_radius, 100)
    
    height = 600
    width = 700
    miny = -0.05
    maxy = target_dist + target_radius + 0.1
    minx = -(width/height) * (maxy - miny) / 2
    maxx = -minx

    scene = GLMakie.lines(
        circle_pts[:, 1], circle_pts[:, 2];
        color=:black, linewidth=3,
        axis=(limits=(minx, maxx, miny, maxy),),
        figure = (resolution = (width, height),)
    )

    longest_traj_idx = argmax(length.(trajs))
    full_wind_traj = trajs[longest_traj_idx].windθs
    num_switches = zeros(Int, length(trajs))

    traj_lines = [GLMakie.Node(GLMakie.Point2f0[]) for _ in 1:length(trajs)]

    windx = GLMakie.Node([0.0])
    windy = GLMakie.Node([-0.3])

    for lines in traj_lines
        GLMakie.lines!(lines, color=:red, linewidth=linewidth, overdraw=true)
    end

    switch_pts = GLMakie.Node(GLMakie.Point2f0[])
    GLMakie.scatter!(switch_pts, color=:blue, markersize=markersize)

    GLMakie.arrows!([-1.0], [1.5], windx, windy, arrowcolor=:green, arrowsize=arrowsize)

    frame_cuts = 1:freq:maximum(length.(trajs))

    GLMakie.record(scene, fname, frame_cuts; framerate=30) do i
        for (t, traj) in enumerate(trajs)
            if i > length(traj)
                continue
            end
            subtraj_len = min(freq, length(traj) - i)
            xs = view(traj.xs, i:i+subtraj_len)
            ys = view(traj.ys, i:i+subtraj_len)
            new_pts = GLMakie.Point2f0.(xs, ys)

            traj_lines[t][] = append!(traj_lines[t][], new_pts)

            # Check if we've hit the next switching point
            if num_switches[t] < length(traj.switch_pts)
                (next_sx, next_sy) = traj.switch_pts[num_switches[t] + 1]
                if any(xs .== next_sx) && any(ys .== next_sy)
                    num_switches[t] += 1
                    switch_pts[] = push!(switch_pts[], GLMakie.Point2f0(next_sx, next_sy))
                end
            end

        end
        windx[] = [0.3 * sin(full_wind_traj[i])]
        windy[] = [-0.3 * cos(full_wind_traj[i])] 
    end
end

function anim_trajs(polar_trajs::Vector{PolarTrajectory}, params::ProblemParams, fname; freq=1, markersize=12, linewidth=4, arrowsize=15)
    trajs = convert.(CartTrajectory, polar_trajs, params.target_dist)
    anim_trajs(trajs, params, fname; freq=freq, markersize=markersize, linewidth=linewidth, arrowsize=arrowsize)
end

function _close(x :: Float64, y :: Float64) :: Bool
    return abs(x - y) < 1e-8
end

const default_starting_pts = [
    (1.8, 0.0, 1),
    (1.8, 0.0, 2),
    (√(1.8^2 + 0.7^2), atan(0.7, 1.7), 1),
    (√(1.8^2 + 0.7^2), atan(0.7, 1.7), 2),
    (√(1.8^2 + 0.7^2), -atan(0.7, 1.7), 1),
    (√(1.8^2 + 0.7^2), -atan(0.7, 1.7), 2),
]