module MesaIO

using CSV
using DataFrames

export read_header, read_history, read_profile, read_profiles_index, MesaRun

"""
Parse a MESA-value token: strip surrounding quotes for strings, else parse as Float64.
"""
function _parse_header_value(tok::AbstractString)
    if startswith(tok, "\"") && endswith(tok, "\"")
        return tok[2:end-1]
    end
    return parse(Float64, tok)
end

"""
    read_header(path) -> NamedTuple

Parse the 3-line scalar header block that starts every MESA history.data /
profile*.data file (column-index row, name row, value row).
"""
function read_header(path::AbstractString)
    lines = open(io -> [readline(io) for _ in 1:3], path)
    names_row = Symbol.(split(lines[2]))
    values_row = _parse_header_value.(split(lines[3]))
    return NamedTuple(zip(names_row, values_row))
end

"""
    read_history(path) -> (header, data)

Read a MESA history.data file. `header` is a NamedTuple of the scalar run
metadata (version, initial params, ...). `data` is a DataFrame with one row
per saved model, columns resolved from the file's own column-name row (never
hardcoded), so it stays valid across different enabled history_columns.list.
"""
function read_history(path::AbstractString)
    header = read_header(path)
    data = CSV.read(path, DataFrame; header=6, delim=' ', ignorerepeated=true, skipto=7)
    return (header=header, data=data)
end

"""
    read_profile(path) -> (header, data)

Read a MESA profile*.data file. Same header/data structure as
[`read_history`](@ref). Note: profile zones are numbered from 1 at the
**surface**, not the center.
"""
function read_profile(path::AbstractString)
    header = read_header(path)
    data = CSV.read(path, DataFrame; header=6, delim=' ', ignorerepeated=true, skipto=7)
    return (header=header, data=data)
end

"""
    read_profiles_index(logs_dir; index_filename="profiles.index") -> DataFrame

Parse `LOGS/profiles.index`, mapping profile numbers to model numbers, so
profile files (named sequentially, not by model number) can be joined
against history.data by model_number.
"""
function read_profiles_index(logs_dir::AbstractString; index_filename::AbstractString="profiles.index")
    path = joinpath(logs_dir, index_filename)
    lines = readlines(path)
    rows = [parse.(Int, split(line)) for line in lines[2:end] if !isempty(strip(line))]
    return DataFrame(
        model_number=[r[1] for r in rows],
        priority=[r[2] for r in rows],
        profile_number=[r[3] for r in rows],
    )
end

"""
    MesaRun(work_dir)

Bundles a MESA work directory with lazily-loaded history/profile data.
Supports holding several runs at once (e.g. for parameter sweeps in
Phase 4) rather than assuming a single hardcoded path.
"""
mutable struct MesaRun
    work_dir::String
    logs_dir::String
    _history::Union{Nothing,NamedTuple}
    _profiles_index::Union{Nothing,DataFrame}

    function MesaRun(work_dir::AbstractString; logs_dir::AbstractString="LOGS")
        return new(work_dir, joinpath(work_dir, logs_dir), nothing, nothing)
    end
end

function history(run::MesaRun)
    if run._history === nothing
        run._history = read_history(joinpath(run.logs_dir, "history.data"))
    end
    return run._history
end

function profiles_index(run::MesaRun)
    if run._profiles_index === nothing
        run._profiles_index = read_profiles_index(run.logs_dir)
    end
    return run._profiles_index
end

"""
    profile(run::MesaRun, profile_number::Integer)

Read `profile<profile_number>.data` from the run's LOGS directory.
"""
function profile(run::MesaRun, profile_number::Integer)
    return read_profile(joinpath(run.logs_dir, "profile$(profile_number).data"))
end

export history, profiles_index, profile

end # module MesaIO
