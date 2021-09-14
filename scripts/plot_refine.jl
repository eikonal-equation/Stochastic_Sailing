using Plots
using Serialization
using SemiLagrangian
using LaTeXStrings
using Debugger

# Make directories to hold all of the outputs, if they don't exist
if !isdir("Plots")
    mkdir("Plots")
end

σlist     = [0.0, 0.05, 0.1, 0.05, 0.05, 0.05]
driftlist = [0.0, 0.0,  0.0, 0.05, 0.1,  0.15]
granularities = [1, 2, 4, 8, 16, 32]

def_grid = initialize_valuegrid(
    SemiLagrangian.default_rbounds,
    SemiLagrangian.default_θbounds,
    SemiLagrangian.default_Δr,
    SemiLagrangian.default_Δθ,
    SemiLagrangian.deterministic_params
)
Nr = def_grid.Nr
Nθ = def_grid.Nθ
Δr = def_grid.Δr
Δθ = def_grid.Δθ

for (σ, drift) in zip(σlist, driftlist)
    # If you have PGFPlotsX.jl installed, you can uncomment out this line
    #  to produce prettier plots.
    # pgfplotsx()

    gs_grids = Array{Float64, 4}(undef, Nθ-1, Nr-1, 2, length(granularities),)
    rwgs_grids = Array{Float64, 4}(undef, Nθ-1, Nr-1, 2, length(granularities),)
    for (i, gran) in enumerate(granularities)
        gs_valgrid = deserialize("Grids/GS_Gran$(gran)Sigma$(σ)Drift$(drift).val")
        rwgs_valgrid = deserialize("Grids/RWGS_Gran$(gran)Sigma$(σ)Drift$(drift).val")
        gs_grids[1:end, 1:end, 1:end, i] .= gs_valgrid.vgrid[1:gran:(Nθ-1)*gran, 1:gran:(Nr-1)*gran, 1:end]
        rwgs_grids[1:end, 1:end, 1:end, i] .= rwgs_valgrid.vgrid[1:gran:(Nθ-1)*gran, 1:gran:(Nr-1)*gran, 1:end]
    end

    gs_grid_diffs = grids .- grids[1:end, 1:end, 1:end, end]
    rwgs_grid_diffs = fastaccgrids .- fastaccgrids[1:end, 1:end, 1:end, end]

    gs_l1_diffs = sum(abs.(gs_grid_diffs), dims=(1, 2, 3)) / ((Nr - 1) * (Nθ - 1) * 2)
    rwgs_l1_diffs = sum(abs.(rwgs_grid_diffs), dims=(1, 2, 3)) / ((Nr - 1) * (Nθ - 1) * 2)
    gs_linf_diffs = maximum(abs.(gs_grid_diffs), dims=(1, 2, 3))
    rwgs_linf_diffs = maximum(abs.(rwgs_grid_diffs), dims=(1, 2, 3))

    gs_l1_diffs     = dropdims(gs_l1_diffs, dims=(1, 2, 3))
    gs_linf_diffs   = dropdims(gs_linf_diffs, dims=(1, 2, 3))
    rwgs_l1_diffs   = dropdims(rwgs_l1_diffs, dims=(1, 2, 3))
    rwgs_linf_diffs = dropdims(rwgs_linf_diffs, dims=(1, 2, 3))
    Δrs = Δr ./ (granularities[end-1:-1:1] ./ 2)
    xticks = Δrs

    plot(Δrs, gs_l1_diffs[end-1:-1:1], 
         marker=true, xlabel=L"$\Delta r$", ylabel=L"$l_1$ Error", label="Gauss-Seidel",
         tickfontsize=14, guidefontsize=14, legendfontsize=14,
         legend=:topleft, xticks=(Δrs, map(string, Δrs)), xscale=:log10, yscale=:log10,)
    plot!(Δrs, rwgs_l1_diffs[end-1:-1:1], marker=true, label="Rowwise Gauss-Seidel")
    p = plot!(Δrs, Δrs, label="", color=:black)
    savefig(p, "Plots/RefineL1CompSigma$(σ)Drift$(drift).png")

    plot(Δrs, gs_linf_diffs[end-1:-1:1],
         marker=true, xlabel=L"$\Delta r$", ylabel=L"$l_\infty$ Error", label="Gauss-Seidel",
         tickfontsize=14, guidefontsize=14, legendfontsize=14,
         legend=:bottomright, xticks=(Δrs, map(string, Δrs)), xscale=:log10, yscale=:log10,)
    plot!(Δrs, rwgs_linf_diffs[end-1:-1:1], marker=true, label="Rowwise Gauss-Seidel")
    p = plot!(Δrs, Δrs ./ Δrs[end], label="", color=:black)
    savefig(p, "Plots/RefineLInfCompSigma$(σ)Drift$(drift).png")

    gr()

    errmap = transpose(abs.(rwgs_grid_diffs[1:end, 1:end, 1, end-1]))
    heatmap(errmap, color=:inferno, proj=:polar, size=(600, 600))
    savefig("Plots/ErrMap_σ$(σ)_D$(drift).png")
end
