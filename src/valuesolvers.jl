function gauss_hermite_switch(ir, iθ, q, valuegrid::ValueGrid, params::ProblemParams) :: Float64
    """ Approximates the expected value function integral at (ir, iθ, q) under wind
         evolution given by the ProblemParams. Uses Gauss-Hermite quadrature with
         n = 3 nodes.
    """
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params
    xs = [-1.22474, 0.0, 1.22474]
    ws = [1/6, 2/3, 1/6] # Absorb the 1/√π factor into these
    qprime = 3 - q
    θ = θbounds[1] + (iθ - 1) * Δθ

    value_estimate = 0.0
    for (x, w) in zip(xs, ws)
        θprime = √(2 * switch_cost) * σ * x + θ + drift * switch_cost
        θprime = mod2pi(θprime)
        value_estimate += w * interpolateθ(ir, θprime, valuegrid, qprime)
    end
    return value_estimate
end

function solve_value(polar::PolarPlot, valuegrid::ValueGrid, params::ProblemParams, λ, ϵ; max_iters=500)
    """The offline phase: perform value iterations to find the value function grid. This 
        version uses an adaptive Δt, but performs the slowest "Gauss-Jacobi"
        iterations where values are not updated in-place during each iteration.
    """
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params

    prev_valuegrid = valuegrid
    cur_valuegrid = deepcopy(valuegrid)
    dirgrid = Array{Float64}(undef, Nθ, Nr, 2)
    switchgrid = zeros(Int8, Nθ, Nr, 2)

    # List of u's to brute force grid search over
    search_us = polar.angles
    
    # At less than this index, we have hit the target
    target_idx = floor(Int, target_radius / Δr) + 1

    # Preallocate an array to hold value function estimates over search space
    Σs = Array{Float64}(undef, size(search_us))
    
    # Value iterations!
    max_diff = Inf
    it = 1
    step_lengths = [1.5 * Δr 1.5 * Δθ]
    @inbounds begin
    while !(max_diff < ϵ)
        p_vgrid = prev_valuegrid.vgrid
        c_vgrid = cur_valuegrid.vgrid

        for q in [1, 2]
            for ir in target_idx+1:Nr
                for iθ in 1:Nθ
                    r = rbounds[1] + (ir - 1) * Δr
                    θ = θbounds[1] + (iθ - 1) * Δθ

                    for i in 1:length(polar)
                        speed_r = polar[i] * dr(r, θ, q, search_us[i])
                        speed_θ = polar[i] * dθ(r, θ, q, search_us[i]) + drift

                        # Have to add a jitter to speeds in case they are 0
                        time_bound_r = 1.5 * Δr / (abs(speed_r) + 1e-8)
                        time_bound_θ = 1.5 * Δθ / (abs(speed_θ) + 1e-8)
                        Δt = min(time_bound_r, time_bound_θ)

                        newr = r + speed_r * Δt
                        newθ = θ + speed_θ * Δt
                        newθ₊ = mod2pi(newθ + √(Δt) * σ)
                        newθ₋ = mod2pi(newθ - √(Δt) * σ)

                        Σs[i] = Δt + exp(-λ * Δt) * 0.5 * (
                            interpolate(newr, newθ₊, prev_valuegrid, q) 
                          + interpolate(newr, newθ₋, prev_valuegrid, q)
                        )
                    end
                    
                    VΔ = gauss_hermite_switch(ir, iθ, q, prev_valuegrid, params) + switch_cost
                    
                    bestΣ_idx = argmin(Σs)
                    bestu = search_us[bestΣ_idx]
                    bestΣ = Σs[bestΣ_idx]
                    dirgrid[iθ, ir, q] = bestu

                    if VΔ < bestΣ
                        max_diff = max(max_diff, abs(p_vgrid[iθ, ir, q] - VΔ))
                        switchgrid[iθ, ir, q] = 1
                        c_vgrid[iθ, ir, q] = VΔ
                    else
                        if !isinf(bestΣ)
                            max_diff = max(max_diff, abs(p_vgrid[iθ, ir, q] - bestΣ))
                        end
                        switchgrid[iθ, ir, q] = 0
                        c_vgrid[iθ, ir, q] = bestΣ
                    end
                end
            end
        end
        
        # Swap prev_valuegrid and cur_valuegrid
        temp = prev_valuegrid
        prev_valuegrid = cur_valuegrid
        cur_valuegrid = temp
        
        max_diff = maximum(filter(!isnan, abs.(cur_valuegrid.vgrid - prev_valuegrid.vgrid)))
        if it % 50 == 0
            println("Iteration $it, Maximum Difference $max_diff")
        end
        if it == max_iters
            println("Stopped after $max_iters iterations")
            return prev_valuegrid
        end
        it += 1
    end
    end

    println("Finished successfully in $it iterations, with final max_diff $max_diff")
    
    # Wrap up grids in the corresponding structs
    dirgrid = DirGrid(dirgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
    switchgrid = SwitchGrid(switchgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
    
    prev_valuegrid, dirgrid, switchgrid, it
end

function solve_value_sweep(polar::PolarPlot, valuegrid::ValueGrid, params::ProblemParams, λ, ϵ; max_iters=500, rowwise=true)
    """The offline phase: perform value iterations to find the value function grid. This 
        version uses an adaptive Δt, and performs the faster "Gauss-Seidel" sweeping
        method. Setting rowwise=true activates the rowwise time discretization.
    """
    valuegrid = deepcopy(valuegrid)
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params

    dirgrid = Array{Float64}(undef, Nθ, Nr, 2)
    switchgrid = zeros(Int8, Nθ, Nr, 2)

    search_us = polar.angles
    
    # Set up boundary values: all distances <= target_radius have exit cost 200
    target_idx = ceil(Int, target_radius / Δr) + 1
    
    # Preallocate an array to hold value function estimates over search space
    Σs = Array{Float64}(undef, length(search_us))

    # Select our discretization strategy
    if rowwise
        Δt_disc = (sr, sθ) -> 1.5 * Δr / (abs(sr) + 1e-8)
    else
        Δt_disc = (sr, sθ) -> min(1.5 * Δr / (abs(sr) + 1e-8), 1.5 * Δθ / (abs(sθ) + 1e-8))
    end

    # Value iterations
    max_diff = Inf
    it = 0
    step_lengths = [1.5 * Δr, 1.5 * Δθ]
    @inbounds begin
    while !(max_diff < ϵ)
        max_diff = 0.0
        for ir in target_idx+1:Nr
            for q in [1, 2]
                for iθ in 1:Nθ
                    r = rbounds[1] + (ir - 1) * Δr
                    θ = θbounds[1] + (iθ - 1) * Δθ

                    for i in 1:length(polar)
                        speed_r = polar[i] * dr(r, θ, q, search_us[i])
                        speed_θ = polar[i] * dθ(r, θ, q, search_us[i]) + drift

                        Δt = Δt_disc(speed_r, speed_θ)

                        newr = r + speed_r * Δt
                        newθ = θ + speed_θ * Δt
                        newθ₊ = mod2pi(newθ + √(Δt) * σ)
                        newθ₋ = mod2pi(newθ - √(Δt) * σ)

                        Σs[i] = Δt + exp(-λ * Δt) * 0.5 * (
                            interpolate(newr, newθ₊, valuegrid, q) 
                          + interpolate(newr, newθ₋, valuegrid, q)
                        )
                    end
                    
                    VΔ = gauss_hermite_switch(ir, iθ, q, valuegrid, params) + switch_cost
                    
                    bestΣ_idx = argmin(Σs)
                    bestu = search_us[bestΣ_idx]
                    bestΣ = Σs[bestΣ_idx]
                    dirgrid[iθ, ir, q] = bestu

                    if VΔ < bestΣ
                        max_diff = max(max_diff, abs(vgrid[iθ, ir, q] - VΔ))
                        switchgrid[iθ, ir, q] = 1
                        vgrid[iθ, ir, q] = VΔ
                    else
                        if !isinf(bestΣ)
                            max_diff = max(max_diff, abs(vgrid[iθ, ir, q] - bestΣ))
                        end
                        switchgrid[iθ, ir, q] = 0
                        vgrid[iθ, ir, q] = bestΣ
                    end
                end
            end
        end
        
        if it % 10 == 0
            println("Iteration $it, Maximum Difference $max_diff")
        end
        if it == max_iters
            println("Stopped after $max_iters iterations")
            dirgrid = DirGrid(dirgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
            switchgrid = SwitchGrid(switchgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
            return valuegrid, dirgrid, switchgrid
        end
        it += 1
    end
    end

    println("Finished successfully in $it iterations, with final max_diff $max_diff")

    # Wrap up grids in the corresponding structs
    dirgrid = DirGrid(dirgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
    switchgrid = SwitchGrid(switchgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
    
    valuegrid, dirgrid, switchgrid, it
end

function solve_value_sweep_nohermite(polar::PolarPlot, valuegrid::ValueGrid, params::ProblemParams, λ, ϵ; max_iters=500, rowwise=true)
    """The offline phase: perform value iterations to find the value function grid. This is
        like `solve_value_sweep`, but uses the inaccurate switching
        operator which does not evolve the wind during the tack switch.
    """
    valuegrid = deepcopy(valuegrid)
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    @unpack target_radius, target_dist, switch_cost, σ, drift, exit_cost = params

    dirgrid = Array{Float64}(undef, Nθ, Nr, 2)
    switchgrid = zeros(Int8, Nθ, Nr, 2)

    search_us = polar.angles
    
    # Set up boundary values: all distances <= target_radius have exit cost 200
    target_idx = ceil(Int, target_radius / Δr) + 1
    
    # Preallocate an array to hold value function estimates over search space
    Σs = Array{Float64}(undef, length(search_us))

    # Select our discretization strategy
    if rowwise
        Δt_disc = (sr, sθ) -> 1.5 * Δr / (abs(sr) + 1e-8)
    else
        Δt_disc = (sr, sθ) -> min(1.5 * Δr / (abs(sr) + 1e-8), 1.5 * Δθ / (abs(sθ) + 1e-8))
    end

    # Value iterations
    max_diff = Inf
    it = 0
    step_lengths = [1.5 * Δr, 1.5 * Δθ]
    @inbounds begin
    while !(max_diff < ϵ)
        max_diff = 0.0
        for ir in target_idx+1:Nr
            for q in [1, 2]
                for iθ in 1:Nθ
                    r = rbounds[1] + (ir - 1) * Δr
                    θ = θbounds[1] + (iθ - 1) * Δθ

                    for i in 1:length(polar)
                        speed_r = polar[i] * dr(r, θ, q, search_us[i])
                        speed_θ = polar[i] * dθ(r, θ, q, search_us[i]) + drift

                        Δt = Δt_disc(speed_r, speed_θ)

                        newr = r + speed_r * Δt
                        newθ = θ + speed_θ * Δt
                        newθ₊ = mod2pi(newθ + √(Δt) * σ)
                        newθ₋ = mod2pi(newθ - √(Δt) * σ)

                        Σs[i] = Δt + exp(-λ * Δt) * 0.5 * (
                            interpolate(newr, newθ₊, valuegrid, q) 
                          + interpolate(newr, newθ₋, valuegrid, q)
                        )
                    end
                    
                    VΔ = vgrid[iθ, ir, 3-q] + switch_cost
                    
                    bestΣ_idx = argmin(Σs)
                    bestu = search_us[bestΣ_idx]
                    bestΣ = Σs[bestΣ_idx]
                    dirgrid[iθ, ir, q] = bestu

                    if VΔ < bestΣ
                        max_diff = max(max_diff, abs(vgrid[iθ, ir, q] - VΔ))
                        switchgrid[iθ, ir, q] = 1
                        vgrid[iθ, ir, q] = VΔ
                    else
                        if !isinf(bestΣ)
                            max_diff = max(max_diff, abs(vgrid[iθ, ir, q] - bestΣ))
                        end
                        switchgrid[iθ, ir, q] = 0
                        vgrid[iθ, ir, q] = bestΣ
                    end
                end
            end
        end
        
        if it % 10 == 0
            println("Iteration $it, Maximum Difference $max_diff")
        end
        if it == max_iters
            println("Stopped after $max_iters iterations")
            dirgrid = DirGrid(dirgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
            switchgrid = SwitchGrid(switchgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
            return valuegrid, dirgrid, switchgrid
        end
        it += 1
    end
    end

    println("Finished successfully in $it iterations, with final max_diff $max_diff")

    # Wrap up grids in the corresponding structs
    dirgrid = DirGrid(dirgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
    switchgrid = SwitchGrid(switchgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
    
    valuegrid, dirgrid, switchgrid, it
end