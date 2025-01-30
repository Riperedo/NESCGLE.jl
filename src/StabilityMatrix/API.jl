include("Grids/grids.jl")
include("StructureFactor/API.jl")

"""
    StabilityMatrix
Input to evaluate the equation of the Self-Consistent Generalized Langevin Equation theory.
# Fields
- `species::Integer`: list of parameters to define the system of interest. The first parameter in the params array must be the volume fraction, and the Second must be the Temperature.
- `params::Array{Float64}`: list of parameters to define the system of interest. The first parameter in the params array must be the volume fraction, and the Second must be the Temperature.
- `S::Function`: Function to generate the static structure factor.
- `k::Array{Float64}`: wave vector grid.
- `system::String`: system identifier.
"""
struct StabilityMatrix
	species::Integer
	params::Array{Float64}
	labels::Vector{String}
	S::Function
	k::Array{Float64}
	system::String
end

# constructors
@doc """Returns the wave vector of an input
"""
wave_vector(sm::StabilityMatrix) = sm.k
@doc """Returns an array with the Static Structure Factor from an input
"""
structure_factor(sm::StabilityMatrix) = sm.S.(sm.k)
@doc """Returns a function to construct the Static Structure Factor from an input
"""
structure_factor_function(sm::StabilityMatrix) = sm.S
@doc """Returns the number of species of the system
"""
number_of_species(sm::StabilityMatrix) = sm.species
@doc """Returns the parameters from an input
"""
parameters(sm::StabilityMatrix) = sm.params
@doc """Returns the volume fraction from an input
"""
volume_fraction(sm::StabilityMatrix) = parameters(sm)[1]
@doc """Returns the kind of system from an input
"""
system_str(sm::StabilityMatrix) = sm.system
@doc """Returns the labels of the system parameters
"""
system_labels(sm::StabilityMatrix) = sm.labels

"""
    SM_HS(ϕ::Float64, k::Array{Float64}; VW = false::Bool)
Constructor of an input for a Hard Sphere system.
# Arguments
- `ϕ::Float64`: volume fraction.
- `k::Array{Float64}`: wave vector array.
# Keywords
- `VW = false::Bool`: Verlet-Weis correction[1].
# References
[1] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
"""
function SM_HS(ϕ::Float64, VW::Float64, k::Array{Float64})
	S = VW == 1.0 ? S_HS_VW(ϕ) : S_HS_PY(ϕ)
	labels = ["phi", "VW"]
	return StabilityMatrix(1, [ϕ, VW], labels, S, k, "HardSphere")
end

function SM_HS(ϕ::Float64, k::Array{Float64})
	return SM_HS(ϕ, 1.0, k)
end

"""
    SM_WCA(ϕ::Float64, T::Float64, k::Array{Float64}; ν = 6::Int64)
Constructor for a soft sphere with pair interaction defined by the Weeks–Chandler–Andersen potential[1].
# Arguments
- `ϕ::Float64`: volume fraction.
- `T::Float64`: Temperature.
- `k::Array{Float64}`: wave vector array.
# Keywords
- `ν = 6::Int64`: Parameter to module the softness of the potential.
# References
[1] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola Phys. Rev. E 87, 052306 – Published 22 May 2013
"""
function SM_WCA(ϕ::Float64, T::Float64, k::Array{Float64}; ν = 6::Int64)
	S = S_WCA_blip(ϕ, T; ν = ν)
	labels = ["phi", "T", "nu"]
	return StabilityMatrix(1, [ϕ, T, ν], labels, S, k, "WCA")
end

"""
    SM_SW(ϕ::Float64, T::Float64, λ::Float64, k::Array{Float64})
Constructor for an input of hard spheres with a Square Well attractive potential.
# Arguments
- `ϕ::Float64`: volume fraction.
- `T::Float64`: Dimensionless temperature.
- `λ::Float64`: range of the well in terms of the diameter sigma.
- `k::Array{Float64}`: wave vector array.
"""
function SM_SW(ϕ::Float64, T::Float64, λ::Float64, k::Array{Float64})
	params = [ϕ, T, λ]
	labels = ["phi", "T", "lambda"]
	S = S_RPA(params; potential = "SquareWell")
	return StabilityMatrix(1, params, labels, S, k, "SquareWell")
end

"""
    SM_Yukawa(ϕ::Float64, ϵ::Float64, z::Float64, k::Array{Float64})
Constructor for an input of hard spheres with a Yukawa attractive potential[1].
# Arguments
- `ϕ::Float64`: volume fraction.
- `ϵ::Float64`: Aplitude potential (effective charge).
- `z::Float64`: Inverse screening length.
- `k::Float64`: Wave vector array.
# References
[1] Yukawa, H. (1935). "On the interaction of elementary particles". Proc. Phys.-Math. Soc. Jpn. 17: 48.
"""
function SM_Yukawa(ϕ::Float64, ϵ::Float64, z::Float64, k::Array{Float64})
	params = [ϕ, ϵ, z]
	labels = ["phi", "epsilon", "z"]
	S = S_RPA(params; potential = "Yukawa")
	return StabilityMatrix(1, params, labels, S, k, "Yukawa")
end

"""
    SM_SALR(ϕ::Float64, ϵ₁::Float64, z₁::Float64, ϵ₂::Float64, z₂::Float64, k::Array{Float64})
Constructor for an input of hard spheres with Sort Attractions and Long Repultions(SALR).
# Arguments
- `ϕ::Float64`: volume fraction.
- `ϵ₁::Float64`: Aplitude potential (effective charge) for the attractive part.
- `z₁::Float64`: Inverse screening length for the attractive part.
- `ϵ₂::Float64`: Aplitude potential (effective charge) for the repulsive part.
- `z₂::Float64`: Inverse screening length for the repulsive part.
- `k::Float64`: Wave vector array.
"""
function SM_SALR(ϕ::Float64, ϵ₁::Float64, z₁::Float64, ϵ₂::Float64, z₂::Float64, k::Array{Float64})
	params = [ϕ, ϵ₁, z₁, ϵ₂, z₂]
	labels = ["phi", "epsilon1", "z1", "epsilon2", "z2"]
	S = S_RPA(params; potential = "SALR")
	return StabilityMatrix(1, params, labels, S, k, "SALR")
end

"""
    SM_AO(ϕ::Float64, ϕₚ⁽ᴿ⁾ ::Float64, ξ::Float64, k::Array{Float64})
Constructor for an input of a mixture of colloid-polimers in the monodisperse approximation [1].
# Arguments
- `ϕ::Float64`: Volume fraction of the colloid.
- `ϕₚ⁽ᴿ⁾::Float64`: Volume fraction of the polimer.
- `ξ::Float64`: Radius of gyration of the polimer.
- `k::Array{Float64}`: wave vector array.
# References
[1] H. N. W. Lekkerkerker et al 1992 EPL 20 559
"""
function SM_AO(ϕ::Float64, ϕₚ⁽ᴿ⁾ ::Float64, ξ::Float64, k::Array{Float64})
	nₚ⁽ᴿ⁾ = 6*ϕₚ⁽ᴿ⁾/(π*(ξ^3))
	params = [ϕ, ϕₚ⁽ᴿ⁾, nₚ⁽ᴿ⁾, ξ]
	labels = ["phiC", "phiP", "np", "xi"]
	S = S_RPA(params, potential = "AsakuraOosawa2")
	return StabilityMatrix(1, params, labels, S, k, "AsakuraOosawa2")
end

"""
    SM_StikyHS(ϕ::Float64, τ::Float64, k::Array{Float64})
Constructor of an input for a Hard Sphere system.
# Arguments
- `ϕ::Float64`: volume fraction.
- `τ::Float64`: .
- `k::Array{Float64}`: wave vector array.
# References
[1]

Contributed by O. Joaquín-Jaime
"""
function SM_StickyHS(ϕ::Float64, τ::Float64, k::Array{Float64})
	S = S_HS_Sticky(ϕ, τ)
	labels = ["phi", "tau"]
	return StabilityMatrix(1, [ϕ, τ], labels, S, k, "StickyHS")
end

"""
    update_SM(sm::StabilityMatrix, params::Array{Float64}) -> StabilityMatrix

Update the stability matrix with new parameters.

# Arguments
- `sm::StabilityMatrix`: An instance of the `StabilityMatrix` type representing the current state of the system.
- `params::Array{Float64}`: The new parameters to update in the stability matrix.

# Returns
- `StabilityMatrix`: A new instance of the `StabilityMatrix` with the updated parameters.

# Example
```julia
sm = StabilityMatrix(
    2, 
    [0.3, 300.0, 1.0, 1.5], 
    ["Volume Fraction", "Temperature", "Param3", "Param4"], 
    x -> 1.0 / (1.0 + x^2), 
    range(0.1, stop=10.0, length=100), 
    "Colloidal Suspension"
)
new_params = [0.4, 320.0, 1.2, 1.6]
updated_sm = update_SM(sm, new_params)
println(updated_sm)
```
"""
function update_SM(sm::StabilityMatrix, params::Array{Float64})
	species = number_of_species(sm)
	k = wave_vector(sm)
	labels = system_labels(sm)
	system = system_str(sm)
	S = structure_factor(params, system)
	return StabilityMatrix(species, params, labels, S, k, system)
end

"""
    update_SM(sm::StabilityMatrix, S::Function) -> StabilityMatrix

Update the structure factor of the stability matrix.

# Arguments
- `sm::StabilityMatrix`: An instance of the `StabilityMatrix` type representing the current state of the system.
- `S::Function`: The new structure factor function to update in the stability matrix.

# Returns
- `StabilityMatrix`: A new instance of the `StabilityMatrix` with the updated structure factor.

# Example
```julia
sm = StabilityMatrix(
    2, 
    [0.3, 300.0, 1.0, 1.5], 
    ["Volume Fraction", "Temperature", "Param3", "Param4"], 
    x -> 1.0 / (1.0 + x^2), 
    range(0.1, stop=10.0, length=100), 
    "Colloidal Suspension"
)
new_S = x -> 1.0 / (1.0 + 2*x^2)
updated_sm = update_SM(sm, new_S)
println(updated_sm)
```
"""
function update_SM(sm::StabilityMatrix, S::Function)
	species = number_of_species(sm)
	params = parameters(sm)
	labels = system_labels(sm)
	k = wave_vector(sm)
	system = system_str(sm)
	return StabilityMatrix(species, params, labels, S, k, system)
end

"""
    Copy(sm::StabilityMatrix) -> StabilityMatrix

Create a copy of the given stability matrix.

# Arguments
- `sm::StabilityMatrix`: An instance of the `StabilityMatrix` type representing the current state of the system.

# Returns
- `StabilityMatrix`: A new instance of the `StabilityMatrix` that is a copy of the given stability matrix.

# Example
```julia
sm = StabilityMatrix(
    2, 
    [0.3, 300.0, 1.0, 1.5], 
    ["Volume Fraction", "Temperature", "Param3", "Param4"], 
    x -> 1.0 / (1.0 + x^2), 
    range(0.1, stop=10.0, length=100), 
    "Colloidal Suspension"
)

sm_copy = Copy(sm)
println(sm_copy)
```
"""
function Copy(sm::StabilityMatrix)
	species = number_of_species(sm)
	params = parameters(sm)
	labels = system_labels(sm)
	S = structure_factor_function(sm)
	k = wave_vector(sm)
	system = system_str(sm)
	return StabilityMatrix(species, params, labels, S, k, system)
end
