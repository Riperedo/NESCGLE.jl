"""`Δη_(k :: Array{Float64}, F :: Array{Float64}, S :: Array{Float64})`
Auxiliar function to compute the shear viscosity relaxation[1][2].
```math
\\Delta\\eta(\\tau; t) = \\frac{k_BT}{60\\pi^2}\\int_0^\\infty dkk^4\\left[\\frac{d}{dk}\\ln S(k;t)\\right]^2\\left[\\frac{F(k,\\tau;t)}{S(k;t)}\\right]^2
```
# Arguments
- `k :: Array{Float64}`: Wave Vector.
- `F :: Array{Float64}`: Intermediate Scarrering function.
- `S :: Array{Float64}`: Statuc structure Factor.
# References
[1] G.Naegele and J. Bergenholtz. Linear viscoelasticity of colloidal mixtures. The Journal of Chemical Physics, 108(23):9893–9904, 1998.
[2] Joaquín-Jaime O., Peredo-Ortiz R. To be published.
"""
function Δη_(k :: Array{Float64}, F :: Array{Float64}, S :: Array{Float64})
	dk = diff(k)
	dS = diff(S)
	dSdk = dS./dk
	integrando = (k[2:end].^4).*(dSdk.^2).*(F[2:end].^2)./(S[2:end].^4)
	return sum(integrando.*dk)/(60*π*π)
end

"""`Δζ_(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, F :: Array{Float64}, Fs :: Array{Float64})`
Friction memory function from the Self-Consistent Generalized Langevin Equation theory[1].
# Arguments
- `ϕ :: Float64`: Volume fraction
- `k :: Array{Float64}`: Wave Vector.
- `S :: Array{Float64}`: Statuc structure Factor.
- `F :: Array{Float64}`: Intermediate Scarrering function.
- `Fs :: Array{Float64}`: Self-Intermediate Scarrering function.
# References
[1] Laura Yeomans-Reyna and Magdaleno Medina-Noyola. Overdamped van hove function of colloidal suspensions. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, pages 3382–3399, September 2000.
"""
function Δζ_(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, F :: Array{Float64}, Fs :: Array{Float64})
	integrando = (k.^4).*((S.-1.0).^2).*F.*Fs./(S.^2)
	dk = diff(k)
	return sum(integrando[2:end].*dk)/(36*π*ϕ)
end

function Σᵢ(Δζ :: Array{Float64}, F :: Matrix{Float64}, n :: Integer)
	Total = zeros(length(F[:,1]))
	for i in 2:(n-1)
		Total = Total + (Δζ[n+1-i]-Δζ[n-i])*F[:, i]
	end
	return Total
end

function ΣD(Δζ :: Array{Float64}, D :: Array{Float64}, n :: Integer)
	Total = 0.0
	for i in 1:n-1
		#Total = Total + (Δζ[n+1-i]-Δζ[n-i])*D[i] + (D[n+1-i]-D[n-i])*Δζ[i]
		Total = Total + D[i]*(Δζ[n-i+1] - Δζ[n-i])
	end
	return Total
end

function ΣΔr2(Δζ :: Array{Float64}, W :: Array{Float64}, n :: Integer)
	msdtsum = 0.0
	for i in 2:n-1
		#msdtsum = msdtsum + (Δζ[i] * (W[n+1-i] - W[n-i]) )
		msdtsum = msdtsum + (Δζ[i-1] * (W[n+1-i] - W[n-i]) )
	end
	return msdtsum
end

"""`F_dump(k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Δζ :: Array{Float64}, dt :: Float64, n :: Integer)`
Auxiliar function to evalaute the evolution of the (Self-)Intermediate scattering function in the Euler integration step.
F_dump(k, n) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
# Arguments
- `k :: Array{Float64}`: Wave Vector.
- `S :: Array{Float64}`: Structure factor.
- `F :: Matrix{Float64}`: (Self-)Intermediate scattering Function.
- `Δζ :: Array{Float64}`: Friction memory function.
- `dt :: Float64`: Time step.
- `n :: Integer`: Current step.
"""
function F_dump(k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Δζ :: Array{Float64}, dt :: Float64, n :: Integer)
	kc = 2*π*1.305
	λ = ((k.^2)/(kc^2) .+ 1).^(-1)
	α = ((k.^2)./S .+ λ*Δζ[1] .+ 1/dt).^(-1)
	return α.*(λ.*(Δζ[n-1]*F[:,1] - Σᵢ(Δζ, F, n)) + F[:,n-1]/dt)
end

"""`F_it(k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Δζ :: Array{Float64}, dt :: Float64, n :: Integer)`
Auxiliar function to evalaute the evolution of the (Self-)Intermediate scattering function in the Euler integration step.
F_it(k, n) = α(k)λ(k)Δζₙ[S(k)-F₁(k)]
# Arguments
- `k :: Array{Float64}`: Wave Vector.
- `S :: Array{Float64}`: Structure factor.
- `F :: Matrix{Float64}`: (Self-)Intermediate scattering Function.
- `Δζ :: Array{Float64}`: Friction memory function.
- `dt :: Float64`: Time step.
- `n :: Integer`: Current step.
"""
function F_it(k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Δζ :: Array{Float64}, dt :: Float64, n :: Integer)
	kc = 2*π*1.305
	λ = ((k.^2)/(kc^2) .+ 1).^(-1)
	α = ((k.^2)./S .+ λ*Δζ[1] .+ 1/dt).^(-1)
	return α.*λ.*(S[:]-F[:,1])*Δζ[n]
end
