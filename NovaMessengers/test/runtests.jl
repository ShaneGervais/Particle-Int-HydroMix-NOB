using Test
using DataFrames
using NovaMessengers

const DATA = joinpath(@__DIR__, "data")

@testset "MesaIO.read_header" begin
    h = read_header(joinpath(DATA, "history_sample.data"))
    @test h.version_number == "r25.12.1"
    @test h.compiler == "gfortran"
    @test h.msun isa Float64
    @test h.msun > 0
end

@testset "MesaIO.read_history" begin
    result = read_history(joinpath(DATA, "history_sample.data"))
    @test result.header.version_number == "r25.12.1"
    df = result.data
    @test nrow(df) > 0
    for col in (:model_number, :star_age, :log_Lnuc, :total_mass_n13,
                :total_mass_o14, :total_mass_o15, :total_mass_f17,
                :total_mass_f18, :total_mass_ne18, :total_mass_ne19)
        @test col in propertynames(df)
    end
    # burst signature: isotope masses should span many orders of magnitude
    # (rise during the burst, negligible before/after)
    @test maximum(df.total_mass_f18) / minimum(df.total_mass_f18) > 1e10
    @test issorted(df.star_age)
end

@testset "MesaIO.read_profile" begin
    result = read_profile(joinpath(DATA, "profile_sample.data"))
    df = result.data
    @test nrow(df) > 0
    @test :zone in propertynames(df)
    @test :mass in propertynames(df)
    # decay (messenger emission) channels
    for r in values(DECAY_REACTIONS)
        @test Symbol("screened_rate_", r) in propertynames(df)
    end
    # formation (synthesis) channels
    for rs in values(FORMATION_REACTIONS), r in rs
        @test Symbol("screened_rate_", r) in propertynames(df)
    end
    @test :raw_rate_r_f17_pg_ne18 in propertynames(df)
    # zones are numbered from 1 at the surface
    @test df.zone[1] == 1
    # screened rates should never be negative
    @test all(df.screened_rate_r_f17_pg_ne18 .>= 0)
    @test all(df.screened_rate_r_n13_wk_c13 .>= 0)
end

@testset "MesaIO.read_profiles_index" begin
    df = read_profiles_index(DATA; index_filename="profiles_index_sample.data")
    @test nrow(df) > 0
    @test propertynames(df) == [:model_number, :priority, :profile_number]
    @test all(df.model_number .> 0)
end

@testset "NuclearDecay.MESSENGER_ISOTOPES" begin
    expected = Set((:n13, :o14, :o15, :f17, :f18, :ne18, :ne19))
    @test Set(keys(MESSENGER_ISOTOPES)) == expected
    @test all(iso.mode == Symbol("beta+") for iso in values(MESSENGER_ISOTOPES))

    # spot-check against ~/mesa-25.12.1/data/rates_data/weak_info.list values
    @test MESSENGER_ISOTOPES[:n13].halflife_s ≈ 597.9
    @test MESSENGER_ISOTOPES[:n13].daughter == :c13
    @test MESSENGER_ISOTOPES[:f18].halflife_s ≈ 7109.281
    @test MESSENGER_ISOTOPES[:ne18].halflife_s ≈ 1.733107

    @test decay_constant(MESSENGER_ISOTOPES[:n13]) ≈ log(2) / 597.9
end

@testset "NuclearDecay.mass_number" begin
    @test mass_number(:f18) == 18
    @test mass_number(:n13) == 13
    @test mass_number(:ne19) == 19
    @test_throws ErrorException mass_number(:notanisotope)
end

@testset "MessengerProduction.DECAY_REACTIONS" begin
    # every messenger isotope has exactly one weak-decay channel -- this IS
    # the messenger emission rate, and must match NuclearDecay's daughter
    for iso in keys(MESSENGER_ISOTOPES)
        @test haskey(DECAY_REACTIONS, iso)
    end
    @test DECAY_REACTIONS[:n13] == :r_n13_wk_c13
    @test DECAY_REACTIONS[:f18] == :r_f18_wk_o18
    @test DECAY_REACTIONS[:ne18] == :r_ne18_wk_f18
end

@testset "MessengerProduction.FORMATION_REACTIONS" begin
    # a deliberately different table from DECAY_REACTIONS -- isotope
    # formation (synthesis), not decay (messenger emission); see the
    # module docs for why the two are not interchangeable
    for iso in keys(MESSENGER_ISOTOPES)
        @test haskey(FORMATION_REACTIONS, iso)
        @test !isempty(FORMATION_REACTIONS[iso])
    end
    # isotopes with multiple formation channels, as found in the CO WD net
    @test length(FORMATION_REACTIONS[:f18]) == 3
    @test length(FORMATION_REACTIONS[:f17]) == 2
    @test length(FORMATION_REACTIONS[:ne19]) == 2
    @test length(FORMATION_REACTIONS[:ne18]) == 1
end

@testset "MessengerProduction._sum_zone_reaction_rates" begin
    import NovaMessengers.MessengerProduction: _sum_zone_reaction_rates

    df = read_profile(joinpath(DATA, "profile_sample.data")).data

    # single-channel isotope (ne18's decay): sum of one column is that column
    @test _sum_zone_reaction_rates(df, (:r_ne18_wk_f18,)) == df.screened_rate_r_ne18_wk_f18

    # multi-channel isotope (f18's formation): sum of three columns
    expected = df.screened_rate_r_o17_pg_f18 .+ df.screened_rate_r_n14_ag_f18 .+ df.screened_rate_r_o15_ap_f18
    @test _sum_zone_reaction_rates(df, FORMATION_REACTIONS[:f18]) ≈ expected

    @test_throws ErrorException _sum_zone_reaction_rates(df, (:r_not_a_real_reaction,))
end

@testset "Transport.escape_probability_neutrino" begin
    @test escape_probability_neutrino() == 1.0
    @test escape_probability_neutrino(1, 2, 3; foo="bar") == 1.0
end

@testset "Transport.klein_nishina_factor" begin
    # Thomson limit: KN factor -> 1 as E -> 0
    @test isapprox(klein_nishina_factor(1e-6), 1.0; atol=1e-4)
    # falls off with energy, always <= 1 in this regime
    @test klein_nishina_factor(0.511) < 1.0
    @test klein_nishina_factor(1.275) < klein_nishina_factor(0.511)
    @test klein_nishina_factor(1.809) < klein_nishina_factor(1.275)
    @test_throws ErrorException klein_nishina_factor(0.0)
    @test_throws ErrorException klein_nishina_factor(-1.0)
end

@testset "Transport.compton_opacity" begin
    # pure helium (X=0) vs pure hydrogen (X=1) at fixed energy: Thomson
    # scaling kappa = 0.2*(1+X) cm^2/g
    E = 0.511
    kn = klein_nishina_factor(E)
    @test compton_opacity(0.0, E) ≈ 0.2 * kn
    @test compton_opacity(1.0, E) ≈ 0.4 * kn
    @test compton_opacity(0.7, E) ≈ 0.2 * 1.7 * kn
end

@testset "Transport.optical_depth_gamma / escape_probability_gamma (synthetic)" begin
    # constant-density shell, radius decreasing linearly from surface (zone 1)
    # to center (zone nz), so tau accumulates by a known amount per zone.
    nz = 101
    r_rsun = collect(range(1.0, 0.0; length=nz))
    rho = fill(1.0, nz)  # g/cm^3
    x_h1 = fill(0.7, nz)
    df = DataFrame(radius=r_rsun, logRho=log10.(rho), x_mass_fraction_H=x_h1)
    rsun_cm = 6.957e10

    tau = optical_depth_gamma(df, 0.511; rsun_cm=rsun_cm)
    @test tau[1] == 0.0
    @test issorted(tau)  # monotonically non-decreasing from surface inward
    @test all(tau[2:end] .> 0)

    # analytic check: uniform rho*kappa over the full shell thickness
    kappa = compton_opacity(0.7, 0.511)
    expected_total_tau = kappa * rho[1] * (r_rsun[1] - r_rsun[end]) * rsun_cm
    @test isapprox(tau[end], expected_total_tau; rtol=1e-6)

    esc = escape_probability_gamma(df, 0.511; rsun_cm=rsun_cm)
    @test esc[1] == 1.0
    @test issorted(esc; rev=true)  # monotonically non-increasing inward
    @test all(0.0 .<= esc .<= 1.0)
    @test esc ≈ exp.(-tau)
end

@testset "Transport.escape_probability_gamma (MesaRun, real data)" begin
    prof = read_profile(joinpath(DATA, "profile_sample.data")).data
    rsun_cm = 6.957e10
    esc = escape_probability_gamma(prof, 0.511; rsun_cm=rsun_cm)
    @test length(esc) == nrow(prof)
    @test all(0.0 .<= esc .<= 1.0)
    @test esc[1] == 1.0  # surface zone always fully transparent to itself
end
