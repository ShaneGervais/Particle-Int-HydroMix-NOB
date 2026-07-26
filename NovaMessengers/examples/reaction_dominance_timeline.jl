using NovaMessengers
using DataFrames
using CairoMakie
using Printf

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
const YEAR_S = 365.25 * 24 * 3600

run = MesaRun(WORK_DIR)

println("Computing energy_breakdown_timeline over all saved profiles (loads each profile once)...")
timeline = energy_breakdown_timeline(run)
age_yr = timeline.age ./ YEAR_S

# --- dominant-reaction handoffs (text) --------------------------------------

println("\n== Dominant reaction handoffs across the burst ==")
prev_dominant = nothing
for i in 1:nrow(timeline)
    row = timeline[i, :]
    ranked = sort([(r, abs(row[r])) for r in TRACKED_REACTIONS]; by=p -> -p[2])
    dominant, rate = ranked[1]
    if dominant !== prev_dominant
        @printf("star_age = %10.4f yr  ->  %-16s  (%.3e erg/s, total = %.3e erg/s)\n",
            age_yr[i], dominant, rate, row.total)
        global prev_dominant = dominant
    end
end

# --- plot: every tracked reaction's |rate| vs star_age, log-y --------------

fig = Figure(size=(900, 550))
ax = Axis(fig[1, 1]; yscale=log10,
    xlabel="Star Age [yr]", ylabel="|Reaction Energy Rate| [erg/s]",
    title="Reaction-resolved energy release through the TNR")

colors = cgrad(:tab20, length(TRACKED_REACTIONS); categorical=true)
for (i, r) in enumerate(TRACKED_REACTIONS)
    rate = abs.(timeline[!, r])
    valid = rate .> 0
    lines!(ax, age_yr[valid], rate[valid]; color=colors[i], label=String(r))
end
lines!(ax, age_yr, timeline.total; color=:black, linewidth=2.5, linestyle=:dash, label="total")

Legend(fig[1, 2], ax; labelsize=8, patchsize=(12, 8), rowgap=1)

outdir = joinpath(@__DIR__, "plt_out")
mkpath(outdir)
outpath = joinpath(outdir, "reaction_dominance_timeline.png")
save(outpath, fig)
println("\nWrote $outpath")
