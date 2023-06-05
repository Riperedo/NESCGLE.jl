""" `βU_SW(T :: Float64, λ :: Float64, k :: Float64)`
Auxiliary function that computes for the Fourier Transform of the perturbation potential.

# Arguments
- `T :: Float64`: Dimensionless temperature.
- `λ :: Float64`: range of the well in terms of the diameter sigma.
- `k :: Float64`: wave vector magnitude.
"""
function βU_SW(T :: Float64, λ :: Float64, k :: Float64)
	k2 = k * k
	if k > 0.0750
		λk = λ * k 
		sink = sin(k)
		cosk = cos(k)
		sinλk = sin(λk)
		cosλk = cos(λk)
		c_aux = ((cosk - λ * cosλk) / k) + ((sinλk - sink) / k2) 
		c_aux = c_aux / k 	
	else 
		λ3 = λ ^ 3
		λ5 = λ3 * λ * λ
		c_aux = (1.0/3.0) * (λ3 - 1.0) - (1.0/30.0) * (λ5 - 1.0) * k2
	end
	return 4.0*π*c_aux/T
end
