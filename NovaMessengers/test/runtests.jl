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
