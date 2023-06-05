""" `decimation(k :: Array{Float64}, S :: Array{Float64}, T :: Array{Float64}, Δζ :: Array{Float64}, D :: Array{Float64}, W :: Array{Float64}, Fs :: Matrix{Float64}, F :: Matrix{Float64}, T_save :: Array{Any}, Δζ_save :: Array{Any}, D_save :: Array{Any}, W_save :: Array{Any}, Fs_save :: Array{Any}, F_save :: Array{Any}, Δη_save :: Array{Any}, index :: Integer)`
Auxiliar function to save and discard the half of the elements of the input arrays to perfrom the caclulus of the next decade in the dynamics code.
# Arguments
- `k :: Array{Float64}`: Wave vector.
- `S :: Array{Float64}`: Structure factor.
- `T :: Array{Float64}`: Waiting Time.
- `Δζ :: Array{Float64}`: Memory friction function.
- `D :: Array{Float64}`: Diffusion coefficient.
- `W :: Array{Float64}`: Mean square displacement.
- `Fs :: Matrix{Float64}`: Self-Intermediate Scattering Function.
- `F :: Matrix{Float64}`: Intermediate Scattering Function.
- `T_save :: Array{Any}`: Save array for waiting time.
- `Δζ_save :: Array{Any}`: Save array for the Memory friction function.
- `D_save :: Array{Any}`: Save array for the Diffusion coefficient.
- `W_save :: Array{Any}`: Save array for the Mean square displacement.
- `Fs_save :: Array{Any}`: Save array for the Self-Intermediate Scattering Function.
- `F_save :: Array{Any}`: Save array for the Intermediate Scattering Function.
- `Δη_save :: Array{Any}`: Save array for shear viscosity relaxation.
- `index :: Integer`: Current Index.
"""
function decimation(k :: Array{Float64}, S :: Array{Float64}, T :: Array{Float64}, 
	Δζ :: Array{Float64}, D :: Array{Float64}, W :: Array{Float64}, Fs :: Matrix{Float64}, 
	F :: Matrix{Float64}, T_save :: Array{Any}, Δζ_save :: Array{Any}, 
	D_save :: Array{Any}, W_save :: Array{Any}, Fs_save :: Array{Any}, 
	F_save :: Array{Any}, Δη_save :: Array{Any}, index :: Integer)
	N = length(T)
	for n in Int(N//2):N
		if n%2 == 1
			append!(T_save, T[n]/2)
			append!(Δζ_save, Δζ[n])
			append!(Fs_save, Fs[index, n])
			append!(F_save, F[index, n])
			append!(Δη_save, Δη_(k, F[:,n], S))
			append!(D_save, D[n])
			append!(W_save, W[n])
		end
	end
	for n in 1:Int(N//2)
		T[n] = T[2*n]
		Δζ[n] = Δζ[2*n]
		D[n] = D[2*n]
		W[n] = W[2*n]
		Fs[:,n] = Fs[:,2*n]
		F[:,n] = F[:,2*n]
	end
	for n in Int(N//2):N
		T[n] = 0.0
		Δζ[n] = 0.0
		D[n] = 0.0
		W[n] = 0.0
		Fs[:,n] .= 0.0
		F[:,n] .= 0.0
	end
	return T, Δζ, D, W, Fs, F, T_save, Δζ_save, D_save, W_save, Fs_save, F_save, Δη_save
end

"""`step_free(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Fs :: Matrix{Float64}, Δζ :: Array{Float64}, T :: Array{Float64}, D :: Array{Float64}, W :: Array{Float64}, dt :: Float64, n :: Integer)`
Funtion that evalautes a free diffusion of a colloidal sytem.
# Arguments
- `ϕ :: Float64`: Volume fraction
- `k :: Array{Float64}`: Wave Vector.
- `S :: Array{Float64}`: Structure factor.
- `F :: Matrix{Float64}`: Intermediate scattering Function.
- `Fs :: Matrix{Float64}`: Self-Intermediate scattering function.
- `Δζ :: Array{Float64}`: Friction memory function.
- `T :: Array{Float64}`: Waiting time.
- `D :: Array{Float64}`: Diffusion coefficient.
- `W :: Array{Float64}` Mean square displacement.
- `dt :: Float64`: Time step.
- `n :: Integer`: Current step.
"""
function step_free(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Fs :: Matrix{Float64}, Δζ :: Array{Float64}, T :: Array{Float64}, D :: Array{Float64}, W :: Array{Float64}, dt :: Float64, n :: Integer)
	α = dt*(k.*k)./S .+ 1.0
	F[:,n+1] = F[:,n]./α
	α = dt*(k.*k) .+ 1.0
	Fs[:,n+1] = Fs[:,n]./α
	Δζ[n+1] = Δζ_(ϕ, k, S, F[:, n+1], Fs[:, n+1])
	T[n+1] = dt*n
	D[n+1] = 1 - Δζ[n+1]*T[n+1]
	W[n+1] = D[n+1]*dt + W[n]
	return T, F, Fs, Δζ, D, W
end

"""`step(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64}, Fs :: Matrix{Float64}, Δζ :: Array{Float64}, T :: Array{Float64}, D :: Array{Float64}, W :: Array{Float64}, dt :: Float64, n :: Integer)`
Funtion that evalautes the colloidal dynamics using an Euler integration step of the form
Fₙ(k) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
		+ α(k)λ(k)Δζₙ[S(k)-F₁(k)]
	 = F_dump(k, n) + F_it(k, n)
where
F_dump(k, n) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
F_it(k, n) = α(k)λ(k)Δζₙ[S(k)-F₁(k)]
and
α(k) = [Δτ⁻¹I + k²DS⁻¹(k) + λ(k)Δζ₁(k)]⁻¹
# Arguments
- `ϕ :: Float64`: Volume fraction
- `k :: Array{Float64}`: Wave Vector.
- `S :: Array{Float64}`: Structure factor.
- `F :: Matrix{Float64}`: Intermediate scattering Function.
- `Fs :: Matrix{Float64}`: Self-Intermediate scattering function.
- `Δζ :: Array{Float64}`: Friction memory function.
- `T :: Array{Float64}`: Waiting time.
- `D :: Array{Float64}`: Diffusion coefficient.
- `W :: Array{Float64}` Mean square displacement.
- `dt :: Float64`: Time step.
- `n :: Integer`: Current step.
"""
function step(ϕ :: Float64, k :: Array{Float64}, S :: Array{Float64}, F :: Matrix{Float64},
	Fs :: Matrix{Float64}, Δζ :: Array{Float64}, T :: Array{Float64}, D :: Array{Float64}, 
	W :: Array{Float64}, dt :: Float64, n :: Integer)
	s = ones(length(k))
	Dump = F_dump(k, S, F, Δζ, dt, n)
	dump = F_dump(k, s, Fs, Δζ, dt, n)
	Δζ[n+1] = Δζ[n]
	F[:,n+1] = Dump + F_it(k, S, F, Δζ, dt, n)
	Fs[:,n+1] = dump + F_it(k, s, Fs, Δζ, dt, n)
	Δζ[n+1] = Δζ_(ϕ, k, S, F[:, n+1], Fs[:, n+1])
	while true
		F[:,n+1] = Dump + F_it(k, S, F, Δζ, dt, n)
		Fs[:,n+1] = dump + F_it(k, s, Fs, Δζ, dt, n)
		Δζ_new = Δζ_(ϕ, k, S, F[:, n+1], Fs[:, n+1])
		if (Δζ[n+1]-Δζ_new)^2 < 1e-10 break end
		Δζ[n+1] = Δζ_new
	end
	T[n+1] = dt*n
	# Computing Self diffusion
	suma = ΣD(Δζ, D, n)
	#d =(1/(1/dt+Δζ[1]))*(-Δζ[Int((n+1 - (n+1)%2)//2)]*D[Int((n+1 - (n+1)%2)//2)] - (Δζ[n+1]-Δζ[n])*D[1] + (1/dt+Δζ[1])*D[n] - suma)
	D[n+1] = (-suma + D[n]/dt)/( (1/dt) + Δζ[1])
	# computing MSD
	#msdtsum = ΣΔr2(Δζ, W, n) - (Δζ[1] + (1/dt)) * W[n]
	W[n+1] = D[n+1]*dt + W[n]
	#W[n+1] = dt*(1 - msdtsum - (Δζ[n+1]*W[1]))/(1 + dt*Δζ[1])
	return T, F, Fs, Δζ, D, W
end

