"""
    PreparationProtocol(process::Vector{Any}, time::Vector{Float64}, str::String, rate::Float64)

Protocol of preparation for the colloidal system.

# Fields
- `process::Vector{Any}` represents the Array of parameters
- `time::Vector{Float64}` represents the Array of time storage
- `str::String` represents the Kind of process
- `rate::Float64` represents the Velocity of the process
"""
struct PreparationProtocol
    process::Vector{Any} # Array of parameters
    time::Vector{Float64} # Array of time storage
    str::String # Kind of process
    rate::Float64 # Velocity of the process
end

"""
    StaticProcess(params::Vector{Float64})

Prepararion protocol to use the Equilibrium version of the NESCGLE Theory.
"""
function StaticProcess(params)
    process = [params]
    time = [0.0]
    return PreparationProtocol(process, time, "StaticProcess", 0.0)
end

"""
    InstantaneousProcess(paramsI::Vector{Float64}, paramsF::Vector{Float64})

Preparation protocol for a sudden change in thermodynamic conditions `paramsI` to `parasmF`.
"""
function InstantaneousProcess(paramsI, paramsF)
    process = [paramsI, paramsF]
    time = [0.0, Inf]
    return PreparationProtocol(process, time, "InstantaneousProcess", Inf)
end

"""
    FiniteRate(paramsI::Vector{Float64}, paramsF::Vector{Float64}, α::Float64, N = 30)

Preparation protocol for a finite rate `α` change in thermodynamic conditions `paramsI` to `parasmF`.
This protocol is an extension of the instantaneous process and is approximated by a sequence of `N` finite instantaneous stepes.
"""
function FiniteRate(paramsI, paramsF, α, N = 30)
    Δ = paramsF.-paramsI
    tp = maximum(abs.(Δ))/α
    Δt = tp/N
    time = [n < N+1 ? n*Δt : Inf for n in 0:N+1]
    process = [n < N+1 ? paramsI + Δ*(n/N) : paramsF for n in 0:N+1]
    return PreparationProtocol(process, time, "FiniteRate", α)
end


"""
    Hysteresis(paramsI::Vector{Float64}, paramsF::Vector{Float64}, α::Float64, N = 30)

Preparation protocol for a hysteresis process at finite rate `α` between the thermodynamic conditions `paramsI` to `parasmF`.
This protocol is an extension of the finite rate process and is approximated by a sequence of `N` finite instantaneous stepes for the initial evolution ramp and another `N` for the final evolutionr ramp.
"""
function Hysteresis(paramsI, paramsF, α, N = 30)
    Δ = paramsF.-paramsI
    tp = maximum(abs.(Δ))/α
    Δt = tp/N
    time_up = [n*Δt for n in 0:N]
    process_up = [paramsI + Δ*(n/N) for n in 0:N]
    time_dn = [n < N+1 ? time_up[end] + n*Δt : Inf for n in 0:N+1]
    process_dn = [n < N+1 ? paramsF - Δ*(n/N) : paramsI for n in 0:N+1]
    time = vcat(time_up, time_dn)
    process = vcat(process_up, process_dn)
    return PreparationProtocol(process, time, "Hysteresis", α)
end


function strip_pp(pp::PreparationProtocol, idx::Integer)
    parameter = []
    for params in pp.process
        append!(parameter, params[idx])
    end
    return parameter
end

function extremal_pp(pp::PreparationProtocol, idx::Integer)
    strip_pp(pp, idx)
    sorted = sort(strip_pp(pp, idx))
    return sorted[1], sorted[end]
end
