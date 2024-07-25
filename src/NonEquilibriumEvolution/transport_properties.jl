"""
	dη(τ::Vector{Any}, Δη::Vector{Any})
Computes the zero shear excess viscosity.
# Attributes
- `τ::Vector{Any}`: Correlation time.
- `Δη::Vector{Any}`: Viscosity shear relaxation.
"""
function dη(τ::Vector{Float64}, ΔG::Vector{Any})
	dη = 0.0
	n = length(τ)
	for i in 3:n
		dlog_τ = log(τ[i]) - log(τ[i-1])
		dη += 0.5*(ΔG[i]*τ[i]+ΔG[i-1]*τ[i-1])*dlog_τ
	end
	return dη
end

"""`b⁻¹(t::Vector{Any}, Δζ::Vector{Any})`
Computes the inverse of the mobility.
# Attributes
- `y::Vector{Any}`: Correlation time.
- `Δζ::Vector{Any}`: Friction memory function.
"""
function b⁻¹(t::Vector{Float64}, Δζ::Vector{Any})
	n = length(t)
	bI = 0.0
	for i in 3:n # t starts at 0.0
		dlog_t = log(t[i]) - log(t[i-1])
		bI += 0.5*(Δζ[i]*t[i]+Δζ[i-1]*t[i-1])*dlog_t
	end
	bI += 1.0
	return bI
end

@doc """Computes the α-relaxation time.
"""
function τα(τ::Vector{Float64}, fs::Vector{Float64})
	idx = fs.>exp(-1)
	if idx == 0
		return τ[1]
	else
		return τ[idx][end]
	end
end
