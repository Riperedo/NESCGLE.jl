# API application programming interface for computing the different structure factors programmed.

# Hard Spheres
include("HardSphere/PercusYevick/PercusYevick.jl")
@doc """`S_HS_PY(ϕ :: Float64)`
Returns a function to construct the Static Structure Factor under the Percus-Yevick Closure[1].
```math
S(k) = \\frac{1}{1 - 24\\phi c(k)}
```
# Arguments
- `ϕ :: Float64`: the volume fraction.
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
function S_HS_PY(ϕ :: Float64)
	f(x) = S_HS_PY(ϕ :: Float64, x)
	return f
end

# Verlet-Weiss correction
include("HardSphere/VerletWeiss/VerletWeiss.jl")
C_HS_VW(ϕ :: Float64, k :: Float64) = C_HS_PY(ϕ_VW(ϕ), k_VW(ϕ, k))
IS_HS_VW(ϕ :: Float64, k :: Float64) = IS_HS_PY(ϕ_VW(ϕ), k_VW(ϕ, k))
S_HS_VW(ϕ :: Float64, k :: Float64) = 1/IS_HS_PY(ϕ_VW(ϕ), k_VW(ϕ, k))
"""`S_HS_VW(ϕ :: Float64)`
Returns a function to compute the Static Structure Factor under the Percus-Yevick Closure[1] using the Verlet-Weis correction[2].
```math
S(k) = \\frac{1}{1 - 24\\phi c(k)}
```
# Arguments
- `ϕ :: Float64`: the volume fraction.
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
[2] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
"""
function S_HS_VW(ϕ :: Float64)
	f(x) = 1/IS_HS_PY(ϕ_VW(ϕ), k_VW(ϕ, x))
	return f
end
"""`S_HS_VW(ϕ :: Float64, k :: Array{Float64})`
Returns aan array that includes the Static Structure Factor under the Percus-Yevick Closure[1] using the Verlet-Weis correction[2].
```math
S(k) = \\frac{1}{1 - 24\\phi c(k)}
```
# Arguments
- `ϕ :: Float64`: the volume fraction.
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
[2] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
"""
function S_HS_VW(ϕ :: Float64, k :: Array{Float64})
	return S_HS_VW(ϕ).(k)
end


# WCA approximation via Blip Function
include("HardSphere/BlipFunction/blip.jl")
function C_WCA_blip(ϕ :: Float64, T :: Float64, k :: Float64; ν = 6)
	λ, λ³ = blip(T, ν=ν)
	return C_HS_VW(λ³*ϕ, λ*k)
end

function IS_WCA_blip(ϕ :: Float64, T :: Float64, k :: Float64; ν = 6)
	λ, λ³ = blip(T, ν=ν)
	return IS_HS_VW(λ³*ϕ, λ*k)
end

"""`S_WCA_blip(ϕ :: Float64, T :: Float64; ν = 6)`
Returns a function to compute the Static Structure Factor under the Percus-Yevick Closure[1] using the Verlet-Weis correction[2] and the blip function[3].
```math
S(k) = \\frac{1}{1 - 24\\phi c(k)}
```
# Arguments
- `ϕ :: Float64`: the volume fraction.
- `T :: Real`: Temperature.
# Keywords
- `ν = 6`: Parameter to module the softness of the potential.
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
[2] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
[3] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola Phys. Rev. E 87, 052306 – Published 22 May 2013
"""
function S_WCA_blip(ϕ :: Float64, T :: Float64; ν = 6)
	λ, λ³ = blip(T, ν=ν)
	f(x) = 1/IS_HS_VW(λ³*ϕ, λ*x)
	return f	
end

"""`S_WCA_blip(ϕ :: Float64, T :: Float64; ν = 6)`
Returns an array that includes the Static Structure Factor under the Percus-Yevick Closure[1] using the Verlet-Weis correction[2] and the blip function[3].
```math
S(k) = \\frac{1}{1 - 24\\phi c(k)}
```
# Arguments
- `ϕ :: Float64`: the volume fraction.
- `T :: Real`: Temperature.
- `k :: Array{Float64}`: wave vector.
# Keywords
- `ν = 6`: Parameter to module the softness of the potential.
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
[2] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
[3] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola Phys. Rev. E 87, 052306 – Published 22 May 2013
"""
function S_WCA_blip(ϕ :: Float64, T :: Float64, k :: Array{Float64}; ν = 6)
	return S_WCA_blip(ϕ, T; ν = 6).(k)
end


# Random Phase Approximation
include("HardSphere/RandomPhaseApproximation/SquareWell/SquareWell.jl")
include("HardSphere/RandomPhaseApproximation/Yukawa/Yukawa.jl")
include("HardSphere/RandomPhaseApproximation/AsakuraOosawa/AsakuraOosawa.jl")

"""`βU(params :: Array{Float64}, k :: Float64, name :: String)`
Auxiliary function to return the Fourier Transform of the pair interaction potential to evaluate the Random Phase Approximation[1].
# Arguments
- `params :: Array{Float64}`: List of parameters to evaluate the potential.
- `k :: Float64`: Wave vector.
- `name :: String`: Potential name selector.
# References
[1] R.V. Sharma and K.C. Sharma. The structure factor and the transport properties of dense fluids having molecules with square well potential are possible generalizations. Physica A: Statistical Mechanics and its Applications, 89(1):213–218, 1977
"""
function βU(params :: Array{Float64}, k :: Float64, name :: String)
	Potentials = Dict(
			"" => 0.0,
			"SquareWell" => βU_SW(params[2], params[3], k),
			"Yukawa" => βU_Yukawa(params[2], params[3], k),
			"SALR" => βU_Yukawa(params[2], params[3], k) - βU_Yukawa(params[4], params[5], k),
			"AsakuraOosawa1" => βU_AO_1(params[1], params[2], params[3], k),
			"AsakuraOosawa2" => βU_AO_2(params[1], params[3], params[4], k) # [ϕ, ϕₚ⁽ᴿ⁾, nₚ⁽ᴿ⁾, ξ]
			)
	return Potentials[name]
end

function C_RPA(params :: Vector{Float64}, k :: Float64; potential = "")
	return C_HS_VW(ϕ, k) + βU(params, k, potential)
end

function IS_RPA(params :: Vector{Float64}, k :: Float64; potential = "")
	ϕ = params[1]
	return 1.0 - (6.0/π)*(4*π*ϕ_VW(ϕ)*C_HS_VW(ϕ, k) + ϕ*βU(params, k, potential))
end

"""`S_RPA(params :: Vector{Float64}, k :: Float64; potential = "")`
Returns the static structure factor using the random phase approximation[1] of the form
```math
S(k) - \\frac{1}{1 - 24\\phi c(k) - \\phi\\beta u(k)}
```
where βu(k) is the Fourier transform of the pair interaction potential.
# Arguments
- `params :: Array{Float64}`: List of parameters to evaluate the potential.
- `k :: Float`: wave vector.
# Keywords
- `potential = "" :: String`: Potential name selector.
# References
[1] R.V. Sharma and K.C. Sharma. The structure factor and the transport properties of dense fluids having molecules with square well potential are possible generalizations. Physica A: Statistical Mechanics and its Applications, 89(1):213–218, 1977
"""
function S_RPA(params :: Vector{Float64}, k :: Float64; potential = "")
	return 1/IS_RPA(params, k, potential = potential)
end

"""`structure_factor(params :: Vector{Float64}, potential = "")`
Returns a function to compute the static structure factor using the random phase approximation[1] of the form
```math
S(k) - \\frac{1}{1 - 24\\phi c(k) - \\phi\\beta u(k)}
```
where βu(k) is the Fourier transform of the pair interaction potential.
# Arguments
- `params :: Array{Float64}`: List of parameters to evaluate the potential.
# Keywords
- `potential = "" :: String`: Potential name selector.
# References
[1] R.V. Sharma and K.C. Sharma. The structure factor and the transport properties of dense fluids having molecules with square well potential are possible generalizations. Physica A: Statistical Mechanics and its Applications, 89(1):213–218, 1977
"""
function S_RPA(params :: Vector{Float64}; potential = "")
	Params = [i<=length(params) ? params[i] : 1.0 for i in 1:5]
	f(x) = 1/IS_RPA(Params, x, potential = potential)
	return f
end

"""`S_RPA(params :: Vector{Float64}, k :: Vector{Float64}; potential = "")`
Returns an array with the static structure factor using the random phase approximation[1] of the form
```math
S(k) - \\frac{1}{1 - 24\\phi c(k) - \\phi\\beta u(k)}
```
where βu(k) is the Fourier transform of the pair interaction potential.
# Arguments
- `params :: Array{Float64}`: List of parameters to evaluate the potential.
- `k :: Vector{Float64}`: Wave vector.
# Keywords
- `potential = "" :: String`: Potential name selector.
# References
[1] R.V. Sharma and K.C. Sharma. The structure factor and the transport properties of dense fluids having molecules with square well potential are possible generalizations. Physica A: Statistical Mechanics and its Applications, 89(1):213–218, 1977
"""
function S_RPA(params :: Vector{Float64}, k :: Vector{Float64}; potential = "")
	return S_RPA(params, potential = potential).(k)
end


include("HardSphere/StickyHS/StickyHS.jl")

@doc """
`S_HS_Sticky(τ::Float64, ϕ::Float64, k::Float64)`
Computes the static structure factor of a Stucky Hard Sphere.
# Arguments
- `τ::Float64`: .
- `ϕ::Float64`: Volume fraction.
- `k::Float64`: Wave vector.
# References
[1] 

Contributed by O. Joquín'Jaime
"""
function S_HS_Sticky(ϕ::Float64, τ::Float64)
	f(x) = S_HS_Sticky(ϕ, τ, x)
	return f
end

"""`structure_factor(params :: Vector{Float64}, system = "")`
Returns a function to construct the static structure factor for a given `system`. The systems that are already programmed are:
* "": HS for PY[1].
* "VW": HS under VW approximation[2].
* "WCA": SS under the blip function approximation[3].
* "SquareWell": Potential under the RPA.
* "Yukawa": Potential under the RPA.
* "SALR": Potential under the RPA.
* "AsakuraOosawa1": Potential under the RPA.
* "AsakuraOosawa2": Potential under the RPA.
* "StickyHS": Sticky Hard Sphere.
# Arguments
- `params :: Vector{Float64}`: list of parameters to evaluate the potential
# Keywords
- `system = ""`: System of interest
# Abbreviations
- HS: Hars Sphere
- PY: Percus-Yevick
- SS: Soft Sphere
- RPA: Random Phase Approximation[4]
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
[2] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
[3] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola Phys. Rev. E 87, 052306 – Published 22 May 2013
[4] R.V. Sharma and K.C. Sharma. The structure factor and the transport properties of dense fluids having molecules with square well potential are possible generalizations. Physica A: Statistical Mechanics and its Applications, 89(1):213–218, 1977
"""
function structure_factor(params :: Vector{Float64}, system = "")
	if system == ""
		return S_HS_PY(params[1])
	elseif system == "VW"
		return S_HS_VW(params[1])
	elseif system == "StickyHS"
		return S_HS_Sticky(params[1], params[2])
	elseif system == "WCA"	
		return S_WCA_blip(params[1], params[2]; ν = params[3])
	elseif system in ["SquareWell", "Yukawa", "SALR", "AsakuraOosawa1", "AsakuraOosawa2"]
		return S_RPA(params, potential = system)
	else
		error("The potential for $potential has not yet been programmed. You must use one of the following:\\n\\t*WCA\\n\\t*VW\\n\\t*SquareWell\\n\\t*Yukawa\\n\\t*SALR\\n\\t*AsakuraOosawa1\\n\\t*AsakuraOosawa2\\nThe default potential returns the Percus Yevick Solution.")
	end
end
