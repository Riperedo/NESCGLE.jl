include("../StabilityMatrix/API.jl")
include("../PreparationProtocol/API.jl")
include("../SelfConsistenScheme/dynamics.jl")
include("../utils/API.jl")
include("transport_properties.jl")
include("utils.jl")

"""
    waiting_times(n::Integer, N::Integer, t0::Float64) -> ExponentialGrid

Generate an exponential grid of waiting times for the NESCGLE computation.

The waiting time is the time after a nonequilibrium process is started, such as a sudden quench in temperature or a slow compression. It is the time after a colloidal system is thermodynamically stressed.

# Arguments
- `n::Integer`: The number of waiting times to generate.
- `N::Integer`: The total number of time steps in the exponential grid.
- `t0::Float64`: The initial waiting time.

# Returns
- `ExponentialGrid`: An exponential grid of waiting times.

# Example
```julia
tw = waiting_times(10, 100, 1.0)
```
"""
function waiting_times(n::Integer, N::Integer, t0::Float64)
    return ExponentialGrid(n, N, t0)
end

waiting_times() = waiting_times(5, 90, 1e-6)


"""
    NESCGLEsolver(sm::StabilityMatrix, pp::PreparationProtocol; kwds...) 

Main function of the NESCGLE system. The two main inputs of this function are the stability matrix and the preparation trotocol.

# Arguments
- `sm::StabilityMatrix` Stability matrix of the problem
- `pp::PreparationProtocol` preparation protocol

# Optional Arguments
- `so=saving_options()` Saving options 
- `tw = waiting_times()` Intermediate waiting times grid
- `path = ""` save path

# Example of usage for a static process
```julia
Nk = 200
kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕ = 0.5
sm = SM_HS(ϕ, k)
pp = StaticProcess(sm.params)
NESCGLEsolver(sm, pp, so = saving_options(true, true, false))
```

# Example of usage for a Instantaneous Process
```julia
Nk = 200
kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕi = 0.5
ϕf = 0.6
sm = SM_HS(ϕi, k)
pp = StaticProcess([ϕi, 1.0], [ϕf, 1.0])
NESCGLEsolver(sm, pp, so = )
```

"""
function NESCGLEsolver(sm::StabilityMatrix, pp::PreparationProtocol; tw = waiting_times(), mute = false)
    @assert length(sm.params) == length(pp.process[1])

    # Dictionary output
    Local = Dict()

    if !mute 
        # Initial evolution
        printstyled("Current process: ", color=:light_green)
        printstyled("$(pp.str)\n", color=:red)
        printstyled("System: ", color=:light_green)
        println("$(sm.system)")
        println("Initial evolution")
    end
    t₀ = 0.0
    u₀ = 0.0
	k = wave_vector(sm)
	Si = structure_factor(sm)
    ϕi = volume_fraction(sm)
	τ, Fs, F, Δζ, ΔG, Δr² = SCGLE(ϕi, k, Si)#; dp = dyn_params(), k_max = 7.2)
    ID = string(0)
    statics = Dict("k"=>k, "S"=>Si)
    dynamics = Dict("tau"=>τ, "Fs"=>Fs, "F"=>F, "DZ"=>Δζ, "DG"=>ΔG, "W"=>Δr²)
    SCS = Dict("Statics"=>statics, "Dynamics"=>dynamics) # Self-Consistent Scheme
    Local[ID] = SCS

    # postprocess
	bI = b⁻¹(τ, Δζ) # mobility
    γI = gammaI(ϕi, k, Si, flag = !mute) # localization length
    Λ, ∂Λ = Lyapunov_Stability(sm, sm)
    if !mute
        printstyled("b⁻¹ = ", color=:light_green)
        printstyled("$bI\n", color=:cyan)
        printstyled("γ⁻¹ = ", color=:light_green)
        printstyled("$γI\n", color=:cyan)
        printstyled("[Λ, ∂Λ] = ", color=:light_green)
        printstyled("$([Λ, ∂Λ/bI])\n", color=:cyan)
    end

    # Saving arrays
    saving_idx = 0
    idx_save = [0]
    t_save = [0.0]
    global_u_save = [0.0]
    local_u_save = [0.0]
    bI_save = [bI]
    η_save = [1 + dη(τ, ΔG)]
    τs_save = [τα(τ, Fs)]
    τc_save = [τα(τ, F)]
    gammaI_save = [γI]
    L_save = [Λ]
    dL_save = [∂Λ/bI]

    # Instantiating target sm
    ℇ = Copy(sm)
    instability = false
    arrest = false
    for idx in 2:length(pp.process)
        # updating StabilityMatrix
        params = pp.process[idx]
        ℇ = update_SM(ℇ, params)

        # Looking for thermodynamic instabilities
		instability = minimum(structure_factor(ℇ)) < 0.0
        if instability && !mute
            printstyled("Warning! ", color=:light_yellow)
            println("Themodynamic Instability")
        end

        # Local initial configuration
        u = 0.0
        t₀ = pp.time[idx-1]
        next_step = false
        for index in eachindex(tw)
            if !mute
                printstyled("\nCurrent process: ", color=:light_green)
                printstyled("$(pp.str)\n", color=:light_red)
                printstyled("System: ", color=:light_green)
                println("$(sm.system)")
                printstyled("Params: ", color=:light_green)
                println("$(sm.labels)")
                printstyled("Process: ", color=:light_green)
                println("$(sm.params) → $(ℇ.params)")
                if instability println("Warning!\nThemodynamic Instability") end
            end
            t = minimum([pp.time[idx-1] + tw[index], pp.time[idx]])
            if !mute
                printstyled("Current time interval: ", color=:light_green)
                println("[$(pp.time[idx-1]), $(pp.time[idx])]")
                printstyled("tw = ", color=:light_green)
                printstyled("$t\n", color=:cyan)
            end
            # Updating time parameters
            Δt = t-t₀
            u₀ += Δt/bI
            u += Δt/bI
            S = NE_structure_factor_function(sm, ℇ, u)
            smU = update_SM(ℇ, S)

            if !mute
                printstyled("u = ", color=:light_green)
                printstyled("$u₀\n", color=:cyan)
            end

            # Local time evolution
            k = wave_vector(smU)
            Su = structure_factor(smU)
            ϕu = volume_fraction(ℇ)
            τ, Fs, F, Δζ, ΔG, Δr² = SCGLE(ϕu, k, Su)#; dp = dyn_params(), k_max = 7.2)

            # postprocess
            bI = b⁻¹(τ, Δζ)
            t₀ = t
            iterations, gammas, system = Asymptotic(ϕu, k, Su, flag = !mute)
            if (system == "Fluid") || (system == "Zero")
                γI = 0.0
            else
                γI = 1/gammas[Int(maximum(iterations))]
                arrest = true
            end

            Λ, ∂Λ = Lyapunov_Stability(smU, ℇ)
            if !mute
                printstyled("b⁻¹ = ", color=:light_green)
                printstyled("$bI\n", color=:cyan)
                printstyled("γ⁻¹ = ", color=:light_green)
                printstyled("$γI\n", color=:cyan)
                printstyled("[Λ, ∂Λ] = ", color=:light_green)
                printstyled("$([Λ, ∂Λ/bI])\n", color=:cyan)    
            end
        
            #saving data
            if pp.time[idx] == Inf
                saving_idx += 1
                ID = string(saving_idx)
                statics = Dict("k"=>k, "S"=>Su)
                dynamics = Dict("tau"=>τ, "Fs"=>Fs, "F"=>F, "DZ"=>Δζ, "DG"=>ΔG, "W"=>Δr²)
                SCS = Dict("Statics"=>statics, "Dynamics"=>dynamics)
                Local[ID] = SCS
                append!(idx_save, saving_idx)
                append!(t_save, t)
                append!(global_u_save, u₀)
                append!(local_u_save, u)
                append!(bI_save, bI)
                append!(η_save, 1 + dη(τ, ΔG))
                append!(τs_save, τα(τ, Fs))
                append!(τc_save, τα(τ, F))
                append!(gammaI_save, γI)
                append!(L_save, Λ)
                append!(dL_save, ∂Λ/bI)
            end # end if

            # Interrupting cases
            if pp.time[idx-1] + tw[index] >= pp.time[idx] 
                next_step = true
            elseif t > 1e12 
                next_step = true
            elseif Fs[end] > exp(-1) 
                next_step = true
            elseif bI > 1e12
                next_step = true
                break
            end # end if

            # next step
            if next_step && (pp.time[idx] != Inf)
                sm = update_SM(smU, S)
                saving_idx += 1
                append!(idx_save, saving_idx)
                append!(t_save, t)
                append!(global_u_save, u₀)
                append!(local_u_save, u)
                append!(bI_save, bI)
                append!(η_save, 1 + dη(τ, ΔG))
                append!(τs_save, τα(τ, Fs))
                append!(τc_save, τα(τ, F))
                append!(gammaI_save, γI)
                append!(L_save, Λ)
                append!(dL_save, ∂Λ/bI)
                break
            end # end if
        end # end for
    end # end for
    if arrest
        if !mute
            println()
            printstyled("Arrest detected, computing asymptotics\n", color=:light_yellow)
        end
        uᵃ = Asymptotic_NE(sm, ℇ, 2*local_u_save[end]; flag= false)
        Sᵃ = NE_structure_factor_function(sm, ℇ, uᵃ)
        smᵃ = update_SM(ℇ, Sᵃ)
        k = wave_vector(smᵃ)
        S = Sᵃ.(k)
        ϕ = volume_fraction(ℇ)
        τ, Fs, F, Δζ, ΔG, Δr² = SCGLE(ϕ, k, S)
        ID = "A"
        statics = Dict("k"=>k, "S"=>S)
        dynamics = Dict("tau"=>τ, "Fs"=>Fs, "F"=>F, "DZ"=>Δζ, "DG"=>ΔG, "W"=>Δr²)
        SCS = Dict("Statics"=>statics, "Dynamics"=>dynamics)
        Local[ID] = SCS
        
        γIᵃ = gammaI(ϕ, k, S; flag = !mute)
        # updating asymptotics
        global_u_save[end] = global_u_save[end] - local_u_save[end] + uᵃ
        local_u_save[end] = uᵃ
        gammaI_save[end] = γIᵃ
        if !mute
            printstyled("uᵃ = ", color=:light_green)
            printstyled("$uᵃ\n", color=:cyan)
            printstyled("γIᵃ = ", color=:light_green)
            printstyled("$γIᵃ\n", color=:cyan)
        end
    end
    Global = Dict("tw"=>t_save, "global_u"=>global_u_save, "local_u"=>local_u_save, "bI"=>bI_save, "eta"=> η_save, "tau_s"=>τs_save, "tau_c"=> τc_save, "gammaI"=> gammaI_save, "L"=> L_save, "dL"=> dL_save, "index"=>idx_save)
    return Dict("local"=>Local, "global"=>Global)
end
