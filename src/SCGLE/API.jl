# API application programming interface for call the function developed for the 
# Self-Consistent Generalized Langevin Equation theory.

include("Grids/grid.jl")
include("StructureFactor/API.jl")
include("Dynamics/API.jl")
include("ArrestDiagram/ArrestLine/Asymptotic.jl")

#############
#	Grids	#
#############

"""`Cheb_plus_U(x₀ :: Float64, x₁ :: Float64, x₂ :: Float64, N₁ :: Int64, N₂ :: Int64)`
Auxiliar function that returns a composed grid by a Chebishev grid and a uniform grid.
# Arguments
- `x₀ :: Float64`: Initial point.
- `x₁ :: Float64`: Middle point.
- `x₂ :: Float64`: Final porint.
- `N₁ :: Int64`: Number of point in the Chebishev grid.
- `N₂ :: Int64`: number of point in the uniform grid.
"""
function Cheb_plus_U(x₀ :: Float64, x₁ :: Float64, x₂ :: Float64, N₁ :: Int64, N₂ :: Int64)
	# first grid
	Chev, w = ChebyshevGrids(x₀, x₁, N₁)
	x = vcat(Chev, collect(x₁:((x₂-x₁)/N₂):x₂))
	return unique!(x)
end

"""`grid_plus_grid(x :: Array{Float64}, y :: Array{Float64})`
Auxiliar fucntion to concat two grids.
# Arguments
- `x :: Array{Float64}`: First grid.
- `y :: Array{Float64}`: Second grid.
"""
function grid_plus_grid(x :: Array{Float64}, y :: Array{Float64}) 
	@assert x[end] <= y[1] "The first grid must have a domain lower than the second one. 
	The last element of the first grid is "*string(x[end])*", and the first element for the second array is "*string(y[1])*"."
	return unique!(vcat(x, y))
end

#############
#	SCGLE	#
#############

"""`Input_SCGLE`
Input to evaluate the equation of the Self-Consistent Generalized Langevin Equation theory.
# Fields
- `params :: Array{Float64}`: list of parameters to define the system of interest. The first parameter in the params array must be the volume fraction, and the Second must be the Temperature.
- `S :: Function`: Function to generate the static structure factor.
- `k :: Array{Float64}`: wave vector grid.
- `system :: String`: system identifier.
"""
struct Input_SCGLE
	params :: Array{Float64}
	S :: Function
	k :: Array{Float64}
	system :: String
end

# constructors
@doc """Returns the wave vector of an input
"""
wave_vector(I :: Input_SCGLE) = I.k
@doc """Returns an array with the Static Structure Factor from an input
"""
structure_factor(I :: Input_SCGLE) = I.S.(I.k)
@doc """Returns a function to construct the Static Structure Factor from an input
"""
structure_factor_function(I :: Input_SCGLE) = I.S
@doc """Returns the parameters from an input
"""
parameters(I :: Input_SCGLE) = I.params
@doc """Returns the volume fraction from an input
"""
volume_fraction(I :: Input_SCGLE) = parameters(I)[1]
@doc """Returns the kind of system from an input
"""
system(I :: Input_SCGLE) = I.system

"""`Input_HS(ϕ :: Float64, k :: Array{Float64}; VW = false :: Bool)`
Constructor of an input for a Hard Sphere system.
# Arguments
- `ϕ :: Float64`: volume fraction.
- `k :: Array{Float64}`: wave vector array.
# Keywords
- `VW = false :: Bool`: Verlet-Weis correction[1].
# References
[1] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
"""
function Input_HS(ϕ :: Float64, k :: Array{Float64}; VW = false :: Bool)
	S = VW ? S_HS_VW(ϕ) : S_HS_PY(ϕ)
	return Input_SCGLE([ϕ], S, k, "HardSphere")
end

"""`Input_WCA(ϕ :: Float64, T :: Float64, k :: Array{Float64}; ν = 6 :: Int64)`
Constructor for a soft sphere with pair interaction defined by the Weeks–Chandler–Andersen potential[1].
# Arguments
- `ϕ :: Float64`: volume fraction.
- `T :: Float64`: Temperature.
- `k :: Array{Float64}`: wave vector array.
# Keywords
- `ν = 6 :: Int64`: Parameter to module the softness of the potential.
# References
[1] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola Phys. Rev. E 87, 052306 – Published 22 May 2013
"""
function Input_WCA(ϕ :: Float64, T :: Float64, k :: Array{Float64}; ν = 6 :: Int64)
	S = S_WCA_blip(ϕ, T; ν = ν)
	return Input_SCGLE([ϕ, T, ν], S, k, "WCA")
end

"""`Input_SW(ϕ :: Float64, T :: Float64, λ :: Float64, k :: Array{Float64})`
Constructor for an input of hard spheres with a Square Well attractive potential.
# Arguments
- `ϕ :: Float64`: volume fraction.
- `T :: Float64`: Dimensionless temperature.
- `λ :: Float64`: range of the well in terms of the diameter sigma.
- `k :: Array{Float64}`: wave vector array.
"""
function Input_SW(ϕ :: Float64, T :: Float64, λ :: Float64, k :: Array{Float64})
	params = [ϕ, T, λ]
	S = S_RPA(params; potential = "SquareWell")
	return Input_SCGLE(params, S, k, "SquareWell")
end

"""`Input_Yukawa(ϕ :: Float64, ϵ :: Float64, z :: Float64, k :: Array{Float64})`
Constructor for an input of hard spheres with a Yukawa attractive potential[1].
# Arguments
- `ϕ :: Float64`: volume fraction.
- `ϵ :: Float64`: Aplitude potential (effective charge).
- `z :: Float64`: Inverse screening length.
- `k :: Float64`: Wave vector array.
# References
[1] Yukawa, H. (1935). "On the interaction of elementary particles". Proc. Phys.-Math. Soc. Jpn. 17: 48.
"""
function Input_Yukawa(ϕ :: Float64, ϵ :: Float64, z :: Float64, k :: Array{Float64})
	params = [ϕ, ϵ, z]
	S = S_RPA(params; potential = "Yukawa")
	return Input_SCGLE(params, S, k, "Yukawa")
end

"""`Input_SALR(ϕ :: Float64, ϵ₁ :: Float64, z₁ :: Float64, ϵ₂ :: Float64, z₂ :: Float64, k :: Array{Float64})`
Constructor for an input of hard spheres with Sort Attractions and Long Repultions(SALR).
# Arguments
- `ϕ :: Float64`: volume fraction.
- `ϵ₁ :: Float64`: Aplitude potential (effective charge) for the attractive part.
- `z₁ :: Float64`: Inverse screening length for the attractive part.
- `ϵ₂ :: Float64`: Aplitude potential (effective charge) for the repulsive part.
- `z₂ :: Float64`: Inverse screening length for the repulsive part.
- `k :: Float64`: Wave vector array.
"""
function Input_SALR(ϕ :: Float64, ϵ₁ :: Float64, z₁ :: Float64, ϵ₂ :: Float64, z₂ :: Float64, k :: Array{Float64})
	params[ϕ, ϵ₁, z₁, ϵ₂, z₂]
	S = S_RPA(params; potential = "SALR")
	return Input_SCGLE(params, S, k, "SALR")
end

"""`Input_AO(ϕ :: Float64, ϕₚ⁽ᴿ⁾  :: Float64, ξ :: Float64, k :: Array{Float64})`
Constructor for an input of a mixture of colloid-polimers in the monodisperse approximation [1].
# Arguments
- `ϕ :: Float64`: Volume fraction of the colloid.
- `ϕₚ⁽ᴿ⁾ :: Float64`: Volume fraction of the polimer.
- `ξ :: Float64`: Radius of gyration of the polimer.
- `k :: Array{Float64}`: wave vector array.
# References
[1] H. N. W. Lekkerkerker et al 1992 EPL 20 559
"""
function Input_AO(ϕ :: Float64, ϕₚ⁽ᴿ⁾  :: Float64, ξ :: Float64, k :: Array{Float64})
	nₚ⁽ᴿ⁾ = 6*ϕₚ⁽ᴿ⁾/(π*(ξ^3))
	params = [ϕ, ϕₚ⁽ᴿ⁾, nₚ⁽ᴿ⁾, ξ]
	S = S_RPA(params, potential = "AsakuraOosawa2")
	return Input_SCGLE(params, S, k, "AsakuraOosawa2")
end

"""`Input_StikyHS(ϕ :: Float64, τ :: Float64, k :: Array{Float64})`
Constructor of an input for a Hard Sphere system.
# Arguments
- `ϕ :: Float64`: volume fraction.
- `τ :: Float64`: .
- `k :: Array{Float64}`: wave vector array.
# References
[1]

Contributed by O. Joaquín-Jaime
"""
function Input_StickyHS(ϕ :: Float64, τ :: Float64, k :: Array{Float64})
	S = S_HS_Sticky(ϕ, τ)
	return Input_SCGLE([ϕ, τ], S, k, "StickyHS")
end


#################
#	Dynamics	#
#################

"""`SCGLE(I :: Input_SCGLE; dt = 1e-10 :: Float64, nT = 6 :: Int64, decimations = 50 :: Int, flag = false)`
Computes the dynamic evolution from a structural input of a colloidal system using the Self-Consistent Generalized Langevin Equation theory[1]. This function returns a set of arrays that contains the correlation time, the self and collective Intermediate Scatering Function, the friction memory function, the shear viscosisty relaxation, the diffusion coefficient and the mean squared displacement.
# Arguments
- `I :: Input_SCGLE`: Structural input.
# Keywords
- `dt = 1e-10 :: Float64`: Initial time step.
- `nT = 6 :: Int64`: Auxiliar number to set the number of intermediate steps.
- `decimations = 50 :: Int`: Number of decades to perform.
- `flag = false`: flag to print internal dynamics procedure.
# References
[1] Laura Yeomans-Reyna and Magdaleno Medina-Noyola. Overdamped van Hove function of colloidal suspensions. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, pages 3382–3399, September 2000.
"""
function SCGLE(I :: Input_SCGLE; dt = 1e-10 :: Float64, nT = 6 :: Int64, 
	decimations = 50 :: Int, flag = false)
	if flag print("Computing...") end
	ϕ = volume_fraction(I)
	k = wave_vector(I)
	S = structure_factor(I)
	τ, F, Fs, Δζ, Δη, D, W = dynamics(ϕ, k, S, dt, nT, decimations)
	if flag println(" Done!") end
	return τ, Fs, F, Δζ, Δη, D, W
end

"""`Asymptotic(I :: Input_SCGLE; flag = true :: Bool)`
Compute the asymptotic dynamics for a given structural input.
# Arguments
- `I :: Input_SCGLE`: Structural input.
# Keywords
- `flag = true :: Bool`: flag to print internal procedure.
"""
function Asymptotic(I :: Input_SCGLE; flag = true :: Bool)
	ϕ = volume_fraction(I)
	k = wave_vector(I)
	S = structure_factor(I)
	return Asymptotic(ϕ, k, S; flag= flag)
end

# TODO DHS, Mixtures
