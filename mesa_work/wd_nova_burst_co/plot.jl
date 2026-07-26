using NovaMessengers
using DataFrames
using DelimitedFiles
using CairoMakie

const WORK_DIR = @__DIR__

"""
    getcol(data, name) -> Vector

Fetch column `name` from a history/profile DataFrame, falling back to the
log10 form MESA stores instead when the plain column isn't enabled
(`log_<name>` in history.data, e.g. `log_Teff`; `log<name>` in profile*.data,
e.g. `logRho`) -- mirrors the convenience mesa_reader gives Python's plot.py.
"""
function getcol(data::AbstractDataFrame, name::AbstractString)
    sym = Symbol(name)
    hasproperty(data, sym) && return data[!, sym]
    for prefix in ("log_", "log")
        logsym = Symbol(prefix * name)
        hasproperty(data, logsym) && return 10.0 .^ data[!, logsym]
    end
    error("column \"$name\" not found (also tried log_$name / log$name)")
end

function read_plot_info(mesa_dir::AbstractString, filename::AbstractString)
    path = joinpath(mesa_dir, "data", "star_data", "plot_info", filename)
    data = readdlm(path)
    return 10.0 .^ data[:, 1], 10.0 .^ data[:, 2]
end

function main()
    mkpath(joinpath(WORK_DIR, "plt_out"))

    run = MesaRun(WORK_DIR)
    _, h = history(run)

    pidx = profiles_index(run)
    last_profile_number = pidx.profile_number[argmax(pidx.model_number)]
    pheader, p = profile(run, last_profile_number)

    mesa_dir = ENV["MESA_DIR"]

    set_theme!(fontsize=10, figure_padding=8)
    fig = Figure(size=(1000, 600))

    gs_left = GridLayout(fig[1, 1])
    gs_right = GridLayout(fig[1, 2])
    colsize!(fig.layout, 1, Relative(1 / 3))
    colsize!(fig.layout, 2, Relative(2 / 3))
    rowgap!(gs_left, 24)
    rowgap!(gs_right, 0)
    colgap!(gs_right, 48)

    ax1 = Axis(gs_left[1, 1]; xscale=log10, yscale=log10,
        xlabel="Effective Temperature [K]", ylabel="Luminosity [L⊙]")
    ax2 = Axis(gs_left[2, 1]; xscale=log10, yscale=log10,
        xlabel="Density [g/cm³]", ylabel="Temperature [K]")

    ax3 = Axis(gs_right[1, 1]; yscale=log10, title="Last 5 Years")
    ax4 = Axis(gs_right[2, 1]; yscale=log10, ylabel="Eff. Temp. [K]")
    ax5 = Axis(gs_right[3, 1]; yscale=log10,
        xlabel="Star Age [yr]", ylabel="Envelope Mass [M⊙]")

    ax6 = Axis(gs_right[1, 2]; xscale=log10, yscale=log10, title="Final Profile")
    ax7 = Axis(gs_right[2, 2]; xscale=log10, yscale=log10,
        ylabel="Specific Power [erg/g/s]")
    ax8 = Axis(gs_right[3, 2]; xscale=log10, yscale=log10,
        xlabel="Exterior Mass Coordinate [M⊙]", ylabel="Density [g/cm³]")

    # --- Panel 1: HR diagram ------------------------------------------------

    Teff = getcol(h, "Teff")
    L = getcol(h, "L")
    lines!(ax1, Teff, L)
    scatter!(ax1, [Teff[end]], [L[end]]; color=:firebrick3, markersize=10)
    ax1.xreversed = true

    # --- Panel 2: T-Rho profile ---------------------------------------------

    Rho = getcol(p, "Rho")
    T = getcol(p, "T")
    lines!(ax2, Rho, T)

    # freeze axis limits from the real T-Rho track (in log space) before
    # adding the wider-ranging reference curves, so they don't rescale the view
    function padded_log_range(values; frac=0.05)
        lo, hi = extrema(values)
        llo, lhi = log10(lo), log10(hi)
        d = (lhi - llo) * frac
        return (10.0^(llo - d), 10.0^(lhi + d))
    end
    xl = padded_log_range(Rho)
    yl = padded_log_range(T)

    psi4_x, psi4_y = read_plot_info(mesa_dir, "psi4.data")
    hburn_x, hburn_y = read_plot_info(mesa_dir, "hydrogen_burn.data")
    heburn_x, heburn_y = read_plot_info(mesa_dir, "helium_burn.data")
    lines!(ax2, psi4_x, psi4_y; color=:lightgrey, linestyle=:dot)
    lines!(ax2, hburn_x, hburn_y; color=:lightgrey, linestyle=:dash)
    lines!(ax2, heburn_x, heburn_y; color=:lightgrey, linestyle=:dash)

    eps_nuc_p = getcol(p, "eps_nuc")
    high_burning = eps_nuc_p .> 1e7
    mid_burning = (eps_nuc_p .> 1e3) .& (eps_nuc_p .< 1e7)
    low_burning = (eps_nuc_p .> 1) .& (eps_nuc_p .< 1e3)
    scatter!(ax2, Rho[high_burning], T[high_burning]; color=:firebrick3,
        markersize=10, label="eps_nuc > 1e7 erg/g/s")
    scatter!(ax2, Rho[mid_burning], T[mid_burning]; color=:darkorange,
        markersize=7, label="eps_nuc > 1e3 erg/g/s")
    scatter!(ax2, Rho[low_burning], T[low_burning]; color=:goldenrod,
        markersize=5, label="eps_nuc > 1 erg/g/s")
    axislegend(ax2; position=:rb, labelsize=8)

    if xl !== nothing
        xlims!(ax2, xl)
        ylims!(ax2, yl)
    end

    # --- Panels 3-5: last 5 years time series -------------------------------

    star_age = getcol(h, "star_age")
    window = 5.0
    mask = star_age .> (maximum(star_age) - window)
    age = star_age[mask]

    LH = getcol(h, "LH")
    lines!(ax3, age, L[mask]; label="L")
    lines!(ax3, age, LH[mask]; linestyle=:dash, label="L_H")
    ax3.ylabel = "Luminosity [L⊙]"
    axislegend(ax3; position=:lt, labelsize=8)
    ax3.xticklabelsvisible = false

    R = getcol(h, "R")
    p4a = lines!(ax4, age, Teff[mask])
    ax4b = Axis(gs_right[2, 1]; yaxisposition=:right, yscale=log10,
        ylabel="Radius [R⊙]")
    linkxaxes!(ax4, ax4b)
    hidespines!(ax4b)
    hidexdecorations!(ax4b)
    p4b = lines!(ax4b, age, R[mask]; linestyle=:dash, color=Cycled(2))
    Legend(gs_right[2, 1], [p4a, p4b], ["T_eff", "R"];
        tellwidth=false, tellheight=false, halign=:left, valign=:top,
        framevisible=false, labelsize=8, margin=(4, 4, 4, 4))
    ax4.xticklabelsvisible = false

    star_mass = getcol(h, "star_mass")
    he_core_mass = getcol(h, "he_core_mass")
    env_mass = star_mass .- he_core_mass
    abs_mdot = getcol(h, "abs_mdot")
    p5a = lines!(ax5, age, env_mass[mask])
    ax5b = Axis(gs_right[3, 1]; yaxisposition=:right, yscale=log10,
        ylabel="|Mdot| [M⊙/yr]")
    linkxaxes!(ax5, ax5b)
    hidespines!(ax5b)
    hidexdecorations!(ax5b)
    p5b = lines!(ax5b, age, abs_mdot[mask]; linestyle=:dash, color=Cycled(2))
    Legend(gs_right[3, 1], [p5a, p5b], ["Delta M_H", "|Mdot|"];
        tellwidth=false, tellheight=false, halign=:left, valign=:center,
        framevisible=false, labelsize=8, margin=(4, 4, 4, 4))

    # --- Panels 6-8: final profile -------------------------------------------

    mass_p = getcol(p, "mass")
    xm = pheader.star_mass .- mass_p

    # ax6/ax7/ax8 are log-x, zoomed to the outer envelope (xlims 1e-10 to
    # 8e-4). xm ~ 0 exactly at the surface zone, which log10 can't represent,
    # and CairoMakie's clip-space interpolation underflows to garbage when a
    # line's data spans many more orders of magnitude than the visible
    # window (here xm reaches ~1.1 Msun at the core) -- so restrict to a
    # window comfortably covering the display range rather than relying on
    # Makie to clip a 10+ decade tail itself.
    valid = 0 .< xm .<= 1e-2
    xm_v = xm[valid]

    iso_labels = [
        ("h1", "¹H"), ("he4", "⁴He"), ("c12", "¹²C"),
        ("n14", "¹⁴N"), ("o16", "¹⁶O"),
    ]
    for (iso, label) in iso_labels
        lines!(ax6, xm_v, getcol(p, iso)[valid]; label=label)
    end
    ylims!(ax6, 4e-5, 1.5)
    axislegend(ax6; position=:lb, labelsize=8)
    ax6.ylabel = "Mass Fraction"
    ax6.xticklabelsvisible = false

    pp_p = getcol(p, "pp")
    cno_p = getcol(p, "cno")
    lines!(ax7, xm_v, eps_nuc_p[valid]; label="eps_nuc")
    lines!(ax7, xm_v, pp_p[valid]; linestyle=:dash, label="eps_pp")
    lines!(ax7, xm_v, cno_p[valid]; linestyle=:dot, label="eps_CNO")
    axislegend(ax7; position=:rt, labelsize=8)
    ylims!(ax7, 1e-6, 1e10)
    ax7.xticklabelsvisible = false

    p8a = lines!(ax8, xm_v, Rho[valid])
    ylims!(ax8, 2e-3, 2e6)
    ax8b = Axis(gs_right[3, 2]; yaxisposition=:right, yscale=log10,
        ylabel="Temperature [K]")
    linkxaxes!(ax8, ax8b)
    hidespines!(ax8b)
    hidexdecorations!(ax8b)
    p8b = lines!(ax8b, xm_v, T[valid]; linestyle=:dash, color=Cycled(2))
    ylims!(ax8b, 2e6, 1.25e8)
    Legend(gs_right[3, 2], [p8a, p8b], ["rho", "T"];
        tellwidth=false, tellheight=false, halign=:right, valign=:top,
        framevisible=false, labelsize=8, margin=(4, 4, 4, 4))

    for ax in (ax6, ax7, ax8)
        ax.xreversed = true
        xlims!(ax, 1e-10, 8e-4)
    end

    outpath = joinpath(WORK_DIR, "plt_out", "wd_nova_burst_grid.svg")
    save(outpath, fig)
    println("Wrote $outpath")
end

main()
