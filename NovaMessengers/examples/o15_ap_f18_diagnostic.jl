using NovaMessengers
using DataFrames
using Printf

const WORK_DIR = joinpath(@__DIR__, "..", "..", "mesa_work", "wd_nova_burst_co")
const YEAR_S = 365.25 * 24 * 3600

run = MesaRun(WORK_DIR)
timeline = energy_breakdown_timeline(run)
age_yr = timeline.age ./ YEAR_S

reaction = :r_o15_ap_f18
rate = timeline[!, reaction]  # signed, NOT abs -- this is the point of this check

println("== $reaction: raw signed rate (erg/s), sign changes, and magnitude vs. total ==\n")
println("min = $(minimum(rate)), max = $(maximum(rate))")
n_negative = count(<(0), rate)
n_positive = count(>(0), rate)
println("negative rows: $n_negative, positive rows: $n_positive, zero rows: $(count(==(0), rate))\n")

# print every row where this reaction is within 3 orders of magnitude of the
# total, or where it's negative -- either is diagnostic, everything else is noise
println(@sprintf("%12s  %14s  %14s  %14s  %10s", "age [yr]", "this reaction", "total", "ratio_to_total", "sign"))
for i in 1:nrow(timeline)
    r = rate[i]
    tot = timeline.total[i]
    ratio = tot == 0 ? NaN : abs(r) / abs(tot)
    if r < 0 || (tot != 0 && ratio > 1e-3)
        @printf("%12.4f  %14.4e  %14.4e  %14.4e  %10s\n",
            age_yr[i], r, tot, ratio, r < 0 ? "NEGATIVE" : "positive")
    end
end
