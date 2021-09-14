function liverun(polar :: PolarPlot, switchgrid :: SwitchGrid, dirgrid :: DirGrid, params :: ProblemParams, starting_pt, Δt; seed=nothing, num_kill=100000)
    """ Computes a (stochastic) live run trajectory from a computed switchgrid.
        switchgrid: The grid of optimal switches produced by a solver
        dirgrid: The grid of optimal directions to steer produced by a solver
        params: Parameters defining the problem
        starting_pt: A tuple (r, θ, q) of the initial condition
        Δt: The timestep using in Eulerian integration
    """
    @unpack sgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = switchgrid
    dgrid = dirgrid.dgrid
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params
    target_idx = floor(Int, (target_radius - rbounds[1]) / Δr) + 1

    if !isnothing(seed)
        seed!(seed)
    end

    r, θ, q = starting_pt
    θ = mod2pi(θ)
    traj = [(r, θ)]
    
    traj_rs = [r]
    traj_θs = [θ]
    wind_traj = [0.0]
    switch_pts = []
    
    nsteps = 1
    # For trajectory recovery, I must track the wind angle and boat angle separately
    windθ = 0
    while r > target_radius
        relθ = mod2pi(θ - windθ)
        ir, iθ = floor(Int, (r - rbounds[1]) / Δr) + 1, floor(Int, (relθ - θbounds[1]) / Δθ) + 1
        
        if r ≥ rbounds[2] - Δr
            println("Too large r, quitting")
            return traj, wind_traj
        end

        # Determine if we should switch tacks.
        # We do it this way: If the cell I'm in has switch = 1 on 3/4 surrounding grid points, I switch
        num_switches = sgrid[iθ, ir, q] + sgrid[iθ+1, ir, q] + sgrid[iθ, ir+1, q] + sgrid[iθ+1, ir+1, q]
        if num_switches ≥ 2
            q = 3 - q
            push!(switch_pts, (r, θ))
            # Hang at this point for an amount of time equal to the switching cost
            num_steps_stopped = floor(Int, switch_cost / Δt)
            append!(traj_rs, fill(r, num_steps_stopped))
            append!(traj_θs, fill(θ, num_steps_stopped))
            # Evolve the wind further this many steps
            for _ in 1:num_steps_stopped
                dW = √Δt * randn()
                windθ += drift * Δt + σ * dW
                push!(wind_traj, windθ)
            end
            nsteps += num_steps_stopped
            continue
        end
        
        # Find direction to steer. Grab the nearest gridpoint, except
        #  close to the target where we need to grab from the next-largest
        #  radial gridpoint as no control exists at the target radius.
        iθ = round(Int, (relθ - θbounds[1]) / Δθ) + 1
        ir = round(Int, (r - rbounds[1]) / Δr) + 1
        ir = max(target_idx+1, ir)
        u = dgrid[iθ, ir, q]

        speed = polar[findfirst(isequal(u), polar.angles)]
        
        # Stochastic Euler integration
        r += speed * dr(r, relθ, q, u) * Δt
        θ += speed * dθ(r, relθ, q, u) * Δt
        dW = √Δt * randn()
        windθ += drift * Δt + σ * dW

        # Wrap angles
        θ = mod2pi(θ)
        windθ = mod2pi(windθ)
        
        push!(traj_rs, r)
        push!(traj_θs, θ)
        push!(wind_traj, windθ)
        
        nsteps += 1
        if nsteps ≥ num_kill
            println("Too long to reach target!")
            return PolarTrajectory(traj_rs, traj_θs, wind_traj, switch_pts)
        end
    end
    
    PolarTrajectory(traj_rs, traj_θs, wind_traj, switch_pts)
end