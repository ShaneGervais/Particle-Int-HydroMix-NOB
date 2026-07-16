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
    @test :screened_rate_r_f17_pg_ne18 in propertynames(df)
    @test :raw_rate_r_f17_pg_ne18 in propertynames(df)
    # zones are numbered from 1 at the surface
    @test df.zone[1] == 1
    # screened rate should never be negative
    @test all(df.screened_rate_r_f17_pg_ne18 .>= 0)
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

@testset "MessengerProduction._finite_diff_production (analytic)" begin
    import NovaMessengers.MessengerProduction: _finite_diff_production

    # pure decay, no production: N(t) = N0*exp(-lambda*t) exactly satisfies
    # dN/dt = -lambda*N, so the reconstructed "production" should be ~0
    # (up to the finite-difference scheme's O((lambda*dt)^2) discretization
    # error -- keep lambda*dt small so that error is negligible).
    lambda = 0.01
    age = collect(0.0:0.01:500.0)
    N = 1.0e10 .* exp.(-lambda .* age)
    r = _finite_diff_production(age, N, lambda)
    @test all(abs.(r.rate) .< 1.0e-4 .* lambda .* maximum(N))

    # constant production P0 into a decaying species reaches the analytic
    # solution N(t) = (P0/lambda)*(1 - exp(-lambda*t)); reconstructed
    # production should recover P0 away from the coarse-step start.
    P0 = 5.0e6
    N2 = (P0 / lambda) .* (1 .- exp.(-lambda .* age))
    r2 = _finite_diff_production(age, N2, lambda)
    @test all(isapprox.(r2.rate[5:end], P0; rtol=1e-2))
end
