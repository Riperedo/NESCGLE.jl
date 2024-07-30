"""
    Lyapunov_Stability(smI::StabilityMatrix, smF::StabilityMatrix) -> Tuple{Float64, Float64}

Calculate the Lyapunov stability and its derivative for a system transitioning from an initial to a final state.

# Arguments
- `smI::StabilityMatrix`: The stability matrix representing the initial state of the system.
- `smF::StabilityMatrix`: The stability matrix representing the final state of the system.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing:
    - `Λ::Float64`: The Lyapunov stability.
    - `dΛ::Float64`: The derivative of the Lyapunov stability.

# Example
```julia
smI = StabilityMatrix(
    2, 
    [0.3, 300.0, 1.0, 1.5], 
    ["Volume Fraction", "Temperature", "Param3", "Param4"], 
    x -> 1.0 / (1.0 + x^2), 
    range(0.1, stop=10.0, length=100), 
    "Colloidal Suspension"
)
smF = StabilityMatrix(
    2, 
    [0.4, 320.0, 1.2, 1.6], 
    ["Volume Fraction", "Temperature", "Param3", "Param4"], 
    x -> 1.0 / (1.0 + 2*x^2), 
    range(0.1, stop=10.0, length=100), 
    "Colloidal Suspension"
)

Λ, dΛ = Lyapunov_Stability(smI, smF)
println("Lyapunov stability: ", Λ)
println("Derivative of Lyapunov stability: ", dΛ)
```
"""
function Lyapunov_Stability(smI::StabilityMatrix, smF::StabilityMatrix)
    k = smI.k
    k² = k.*k
    k⁴ = k².*k²
    Δk = k[2]-k[1]
    prefactor = Δk/(2*π*π)
    σ = structure_factor(smI)
    ℇ = structure_factor(smF).^(-1)
    Λ = -0.5*prefactor*sum(k².*(σ.*ℇ .- 1.0))
    dΛ = prefactor*sum(k⁴.*(σ.*ℇ .- 1.0).*ℇ)
    return Λ, dΛ
end

"""
	NE_structure_factor(smI::StabilityMatrix, smF::StabilityMatrix, u::Float64)
Returns a function to evaluate the Non-Equilibrium Static structure factor[1]
# Arguments
- `smI::StabilityMatrix`: Initial Structural input.
- `smF::StabilityMatrix`: Final Structural input.
- `u::Float64`: Inner clock.
# References
[1] Pedro Ramírez-González and Magdaleno Medina-Noyola. General nonequilibrium theory of colloid dynamics. Phys. Rev. E, 82:061503, Dec 2010.
"""
function NE_structure_factor_function(smI::StabilityMatrix, smF::StabilityMatrix, u::Float64)
	Si = structure_factor_function(smI)
	Sf = structure_factor_function(smF)
	f(x) = Si(x)*exp(-2*u*(x^2)/Sf(x)) + Sf(x)*(1-exp(-2*u*(x^2)/Sf(x)))
	return f
end


function NE_structure_factor(smI::StabilityMatrix, smF::StabilityMatrix, u::Float64, k::Vector{Float64})
	return NE_structure_factor_function(smI, smF, u).(k)
end

function complete_str(x)
    str = string(x)
    while length(str) < 8
        str *= "0"
    end
    return str[1:8]
end

"""
	Asymptotic_NE(smI::StabilityMatrix, smF::StabilityMatrix, u::Float64; flag= false)
Compute the asymptotic dynaimics for a given non-equilibrium structural input. This function finds an ergodic state.
# Arguments
- `smI::StabilityMatrix`: Initial Structural input.
- `smF::StabilityMatrix`: Final Structural input.
- `u::Float64`: Inner clock.
# Keywords
- `flag = false`: flag to print internal procedure.
"""
function Asymptotic_NE(smI::StabilityMatrix, smF::StabilityMatrix, u::Float64; flag= false)
	ϕ = volume_fraction(smF)
	function condition(U)
		if flag println("Computing u = $U") end
		k = wave_vector(smI)
        S = NE_structure_factor_function(smI, smF, U).(k)
		system = "Dump"
		if minimum(S) <= 0.0
			system = "Zero"
		else
			iterations, gammas, system = Asymptotic(ϕ, k, S, flag = flag)
		end
		return system != "Fluid"
	end
	return bisection(condition, u, 0.0, 1e-10, flag = false)
end

"""
    make_directories(sm::StabilityMatrix, pp::PreparationProtocol) -> String

Create a nested directory structure based on the stability matrix and preparation protocol parameters.

# Arguments
- `sm::StabilityMatrix`: An instance of the `StabilityMatrix` type containing the system's stability parameters.
- `pp::PreparationProtocol`: An instance of the `PreparationProtocol` type containing the preparation protocol details.

# Returns
- `String`: The path to the final directory created.

# Example
```julia
sm = StabilityMatrix(2, [0.3, 300.0, 1.0, 1.5], ["Volume Fraction", "Temperature", "Param3", "Param4"], x -> 1.0 / (1.0 + x^2), range(0.1, stop=10.0, length=100), "Colloidal Suspension")
pp = PreparationProtocol("Protocol1", 0.5, 1.0, 100)

final_path = make_directories(sm, pp)
println("Final path: ", final_path)
```
"""
function make_directories(sm::StabilityMatrix, pp::PreparationProtocol)
    path = make_directory(pp.str) # preparation protocol\
    path = make_directory(path*sm.system) # system\
    if pp.rate > 0.0 && pp.rate < Inf
        path = make_directory(path*"rate_"*complete_str(pp.rate))
    end
    for idx in eachindex(sm.params)
        label = sm.labels[idx]
        if label in exception_dir continue end
        param_1, param_end = extremal_pp(pp,idx)
        if param_1 == param_end
            path = make_directory(path*label*"_"*complete_str(param_end))
        else
            if sm.params[idx] == param_1
                path = make_directory(path*label*"_"*complete_str(param_1)*"_to_"*complete_str(param_end))
            else
                path = make_directory(path*label*"_"*complete_str(param_2)*"_to_"*complete_str(param_1))
            end
        end
    end
    return path
end

@doc"""
    saving_options(SF::Bool, TP::Bool, tw::Bool, folder_creation::Bool)

A type representing the options for saving different aspects of the NESCGLE computation results.

# Fields
- `SF::Bool`: Indicates whether to save the structure factor.
- `TP::Bool`: Indicates whether to save the transport properties.
- `tw::Bool`: Indicates whether to save the waiting times.
- `folder_creation::Bool`: Indicates whether to create a folder for saving the results.

# Example
```julia
options = saving_options(true, false, true, true)
```
"""
struct saving_options 
    SF::Bool # Structure factor
    TP::Bool # Transport properties
    tw::Bool # Waiting times
    folder_creation::Bool # folder creation
end

saving_options() = saving_options(true,true,true,true)

function saving_options(SF::Bool, TP::Bool, tw::Bool)
    folder_creation = SF || TP || tw
    return saving_options(SF, TP, tw, folder_creation)
end


function save_files(sol, sm::StabilityMatrix, pp::PreparationProtocol; so=saving_options(), path = "")
    # making saving folders
    if so.folder_creation && path == ""
        path = make_directories(sm, pp)
    elseif so.folder_creation && path != ""
        make_directory(path)
    end # end if

    for key in keys(sol["local"])
        #keys(sol["local"][key]) : ["Dynamics", "Statics"]
        #keys(sol["local"][key]["Statics"]) : ["S", "k"]
        #keys(sol["local"][key]["Dynamics"]) : ["DZ", "W", "DG", "Fs", "tau", "F"]
        if so.SF 
            statics = sol["local"][key]["Statics"]
            k = statics["k"]
            S = statics["S"]
            save_data(path*"SF_"*key*".dat", [k S], header = "k\tS", flag = false)
        end
        if so.TP 
            dynamics = sol["local"][key]["Dynamics"]
            τ = dynamics["tau"]
            Fs = dynamics["Fs"]
            F = dynamics["F"]
            Δζ = dynamics["DZ"]
            ΔG = dynamics["DG"]
            Δr² = dynamics["W"]
            save_data(path*"TP_"*key*".dat", [τ Fs F Δζ ΔG Δr²], header = "τ\tFs\tF\tΔζ\tΔG\tΔr²", flag = false)
        end
    end

    if so.tw
        #keys(sol["global"]):["local_u", "dL", "L", "tw", "tau_s", "tau_c", "global_u", "eta", "index", "gammaI", "bI"]
        t = sol["global"]["tw"]
        global_u = sol["global"]["local_u"]
        local_u = sol["global"]["global_u"]
        bI = sol["global"]["bI"]
        η = sol["global"]["eta"]
        τs = sol["global"]["tau_s"]
        τc = sol["global"]["tau_c"]
        γI = sol["global"]["gammaI"]
        L = sol["global"]["L"]
        dL = sol["global"]["dL"]
        idx = sol["global"]["index"]
        save_data(path*"waiting_times.dat", [t global_u local_u bI η τs τc γI L dL idx], header = "t\tglobal_u\tlocal_u\tbI\tη\tτs\tτc\tγI\tΛ\t∂Λ\tidx", flag = false)
    end
end
