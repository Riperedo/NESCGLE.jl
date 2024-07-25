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
