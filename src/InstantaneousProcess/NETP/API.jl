# NETP stands for Non-Equilibirum transport properties

"""`dη(τ :: Array{Any}, Δη :: Array{Any})`
Computes the zero shear excess viscosity.
# Attributes
- `τ :: Array{Any}`: Correlation time.
- `Δη :: Array{Any}`: Viscosity shear relaxation.
"""
function dη(τ :: Array{Any}, Δη :: Array{Any})
	dη = 0.0
	n = length(τ)
	for i in 2:n
		dlog_τ = log(τ[i]) - log(τ[i-1])
		dη += 0.5*(Δη[i]*τ[i]+Δη[i-1]*τ[i-1])*dlog_τ
	end
	return dη
end

"""`b⁻¹(t :: Array{Any}, Δζ :: Array{Any})`
Computes the inverse of the mobility.
# Attributes
- `y :: Array{Any}`: Correlation time.
- `Δζ :: Array{Any}`: Friction memory function.
"""
function b⁻¹(t :: Array{Any}, Δζ :: Array{Any})
	n = length(t)
	bI = 0.0
	for i in 2:n
		dlog_t = log(t[i]) - log(t[i-1])
		bI += 0.5*(Δζ[i]*t[i]+Δζ[i-1]*t[i-1])*dlog_t
	end
	bI += 1.0
	return bI
end

@doc """Computes the α-relazation time.
"""
τα(τ :: Array{Any}, fs :: Array{Any}) = τ[fs.>exp(-1)][end]
