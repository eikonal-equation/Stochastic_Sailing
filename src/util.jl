@inline function interpolate(r, θ, valuegrid::ValueGrid, q) :: Float64
    """ Computes a bilinear interpolation at (r, θ) of the values of the four nearest gridpoints.
    """
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    grid = @view vgrid[:, :, q]

    # 1e-8 jitter is to handle case when r is exactly a multiple of Δr, won't affect others
    ir = floor(Int, (r - rbounds[1]) / Δr + 1e-14) + 1
    iθ = floor(Int, (θ - θbounds[1]) / Δθ + 1e-14) + 1
    frac_r = (r - rbounds[1] - Δr * (ir - 1)) / Δr
    frac_θ = (θ - rbounds[1] - Δθ * (iθ - 1)) / Δθ
    
    # Add some jitter to prevent 0 * Inf. Not needed if not using Infs
    frac_r = clamp(frac_r, 1e-14, 1.0 - 1e-14)
    frac_θ = clamp(frac_θ, 1e-14, 1.0 - 1e-14)
    
    # We have overshot the target: Just return boundary value
    if ir ≤ 0
        return grid[1, 1]
    end
    # We have gone backwards off the grid. Heuristic: Return Inf
    if ir >= Nr
        return Inf
    end
    
    val_bl = grid[iθ, ir]
    val_br = grid[iθ, ir+1]
    val_tl = grid[iθ+1, ir]
    val_tr = grid[iθ+1, ir+1]

    # Linearly interpolate along top and bottom edges, then along r axis between new interpolated points
    # This form should have better numerical stability than p1 + (p2 - p1) * frac, but needs 2 multiplies
    interp_top = frac_r * val_tr + (1.0 - frac_r) * val_tl
    interp_bot = frac_r * val_br + (1.0 - frac_r) * val_bl
    return frac_θ * interp_top + (1.0 - frac_θ) * interp_bot
end

@inline function interpolateθ(ir, θ, valuegrid::ValueGrid, q) :: Float64
    """ Like interpolate, but for the special case of gauss_hermite_switch where I know
          that r falls on a grid value
    """
    @unpack vgrid, Nr, Nθ, rbounds, θbounds, Δr, Δθ = valuegrid
    grid = @view vgrid[:, :, q]

    iθ = floor(Int, (θ - θbounds[1]) / Δθ + 1e-14) + 1
    frac_θ = (θ - rbounds[1] - Δθ * (iθ - 1)) / Δθ

    return (1.0 - frac_θ) * grid[iθ, ir] + frac_θ * grid[iθ+1, ir]
end

# Initializes a valuegrid with the boundary conditions defined by params
function initialize_valuegrid(rbounds, θbounds, Δr, Δθ, params::ProblemParams) :: ValueGrid
    Nr = ceil(Int, (rbounds[2] - rbounds[1]) / Δr) + 1
    Nθ = ceil(Int, (θbounds[2] - θbounds[1]) / Δθ) + 1
    grid = Array{Float64}(undef, Nθ, Nr, 2)
    
    target_idx = floor(Int, (params.target_radius - rbounds[1]) / Δr) + 1
    grid[:, 1:target_idx, :] .= params.exit_cost
    grid[:, target_idx+1:end, :] .= 1e6

    ValueGrid(grid, Nr, Nθ, rbounds, θbounds, Δr, Δθ)
end

const deterministic_params = ProblemParams(0.1, 1.8, 2.0, 0.0, 0.0, 0.0)
const stochastic_params = ProblemParams(0.1, 1.8, 2.0, 0.1, 0.0, 0.0)
const drift_stochastic_params = ProblemParams(0.1, 1.8, 2.0, 0.1, 0.1, 0.0)
const default_rbounds = (0.0, 2.0)
const default_θbounds = (0.0, 2π)
const default_Δr = 0.02
const default_Δθ = 0.05

default_grid() = initialize_valuegrid(
    default_rbounds,
    default_θbounds,
    default_Δr,
    default_Δθ,
    deterministic_params
)

finer_grid() = initialize_valuegrid(
    default_rbounds,
    default_θbounds,
    default_Δr / 2,
    default_Δθ / 2,
    deterministic_params
)

extrafine_grid() = initialize_valuegrid(
    default_rbounds,
    default_θbounds,
    default_Δr / 4,
    default_Δθ / 4,
    deterministic_params
)

thefinestgrid() = initialize_valuegrid(
    default_rbounds,
    default_θbounds,
    default_Δr / 8,
    default_Δθ / 8,
    deterministic_params 
)

# These are the time derivatives in the differential equation, missing a factor of the boat speed
# r - dist, θ - angle, q - tack, u - steering angle
@inline dr(r, θ, q, u) = -cos(θ - (-1)^q * u)
@inline dθ(r, θ, q, u) = sin(θ - (-1)^q * u) / r

# Versions where (r, θ) are bundled together
@inline dr((r, θ), q, u) = dr(r, θ, q, u)
@inline dθ((r, θ), q, u) = dθ(r, θ, q, u)