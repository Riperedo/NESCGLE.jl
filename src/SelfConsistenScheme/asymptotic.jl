"""`γ⁻¹(ϕ::Float64, k::Array{Float64}, S::Array{Float64}, γ::Float64)`
Fucntion to compute the inverse of the localization length.
# Arguments
- `ϕ::Float64`: Volume fraction.
- `k::Array{Float64}`: wave vector.
- `S::Array{Float64}`: Structure factor.
- `γ::Float64`: Seed value.
"""
function γ⁻¹(ϕ::Float64, k::Array{Float64}, S::Array{Float64}, γ::Float64)
	kc = 1.305*2*π
	λ = ((k/kc).^2 .+ 1.0).^(-1)
	fs = λ./(λ .+ (k.^2)*γ)
	f = λ.*S./(λ.*S .+ (k.^2)*γ)
	integrand = (k.^4).*(S.-1).*f.*fs.*(S.-1)./S
	dk = diff(k)
	I = sum(integrand[2:end].*dk)
	return max(I/(36.0*π*ϕ), 1e-12)
end

function selector(γ::Float64, M::Float64, Error_Relativo::Float64, infinito::Float64)
	Converge = abs(M) < Error_Relativo 
	Diverge = γ > infinito

	if Diverge
		return "Fluid"
	elseif Converge
		return "Glass"
	else
		return "Dump"
	end
end

"""`Asymptotic(ϕ, k::Array{Float64}, S::Array{Float64}; flag= false)`
Function to determine the kind of system studed. The possible outputs can be:
*"Dump": for no-confident result.
* "Fluid": For a ergodic state.
* "Glass": For an arrested state.
* "Zero": For thermodynamic instabilities.
# Arguments
- `ϕ::Float64`: Volume fraction.
- `k::Array{Float64}`: wave vector.
- `S::Array{Float64}`: Structure factor.
# Keywords
- `flag = false`: Option to print the internal computing steps.
"""
function Asymptotic(ϕ, k::Array{Float64}, S::Array{Float64}; flag= false)
	# itertions
	It_MAX = 1000
	decimo = div(It_MAX, 50)
	if flag
		printstyled("|-------------------------|------------------------| <- 100%\n", color=:green)
		printstyled("|", color=:green)
		#println("|-------------------------|------------------------| <- 100%")
		#print("|")
	end

	#outputs
	system = "Dump"
	iterations = zeros(It_MAX+1)
	gammas = zeros(It_MAX+1)

	# seed
	γ = 0.00001
	gammas[1] = γ

	# Main loop
	it = 1
	while true
		if it % decimo == 0 && flag
		    printstyled("#", color=:green)
		end
		#system= Dump()
		γ_new = γ⁻¹(ϕ, k, S, γ)
		if γ_new > 0.0
			γ_new = 1/γ_new
			convergencia = (γ - γ_new)/γ_new
		    γ = γ_new
		else
			printstyled(γ_new, color=:red)
			system = "Zero"
			break
		end
		
		iterations[it+1] = it
	
		gammas[it+1] = γ
		
		system = selector(γ, convergencia, 0.0001, 1e10)
		if system == "Fluid"
			if flag printstyled(system, color=:light_blue) end
			break
        elseif system == "Glass"
			if flag printstyled(system, color=:light_red) end
			break
		end
		it += 1
		if it > It_MAX break end
	end
	if it < It_MAX && flag
		qwerty = div(It_MAX-it, decimo)
		for dot in 1:qwerty-4 print(" ") end
	end
	if flag printstyled("| Done!\n", color=:green) end
	return iterations, gammas, system
end

"""`Asymptotic_structure(ϕ::Float64, k::Array{Float64}, S::Array{Float64})`
Function to compute the asymptotic static structure factor.
# Arguments
- `ϕ::Float64`: Volume fraction.
- `k::Array{Float64}`: wave vector.
- `S::Array{Float64}`: Structure factor.
"""
function Asymptotic_structure(ϕ::Float64, k::Array{Float64}, S::Array{Float64})
	iterations, gammas, system = Asymptotic(ϕ, k, S)
	γ = gammas[Int(maximum(iterations))]
	kc = 1.305*2*π
	λ = ((k/kc).^2 .+ 1.0).^(-1)
	F∞ = λ.*S.*S./(λ.*S .+ (k.^2)*γ)
	Fs∞ = λ./(λ .+ (k.^2)*γ)
	integrando = (k.^4).*((S.-1).^2).*F∞.*Fs∞./(S.^2)
	dk = diff(k)
	Δζ∞ = sum(integrando.*dk)/(36.0*π*ϕ)
	return Δζ∞, Fs∞, F∞
end

"""`bisection(condition::Function, A::Real, T::Real, tolerance::Real; flag = false)`
Function to perform a sucessive bisection between two interest points.

According to Zeno's paradox, Achilles can never catch up with the tortoise because he has an infinite number of finite catch-ups to make. This apparent paradox is used to argue that motion is impossible and is simply an illusion. When Achilles reaches the point where the tortoise started the race, the tortoise will have already moved on. The paradox concerns a race between the fleet-footed Achilles and a slow-moving tortoise. The two start moving at the same moment, but if the tortoise is initially given a head start and continues to move ahead, Achilles can run at any speed and will never catch up with it.

# Arguments
- `condition::Function`: Criteria to perform a step.
- `A::Real`: Initial condition.
- `T::Real`: Maximum position of the target.
- `tolerance::Real`: Tolerance. 
# Keywords
- `flag = false`: flag to print the internal computing of the procedure.
"""
function bisection(condition::Function, A::Real, T::Real, tolerance::Real; 
	flag = false)
	Achilles = A
	Tortoise = T
	δ = 0.5
	paso(t, a, Δ) = Δ*(t - a)
	while Achilles != Tortoise
		δA = paso(Tortoise, Achilles, δ)
		if condition(Achilles+δA) # Si la condition se cumple se acepta el paso
			Achilles += δA
		else
			δ *= 0.5
		end

		if abs(δA) < tolerance
			break
		end
		if flag println(Achilles, " ", δA, " ", Achilles+δA) end
	end
	return Achilles
end

function gammaI(ϕ::Float64, k::Vector{Float64}, S::Vector{Float64}; flag = false)
    iterations, gammas, system = Asymptotic(ϕ, k, S, flag = flag)
	if (system == "Fluid") || (system == "Zero")
		return 0.0
	else
		return 1/gammas[Int(maximum(iterations))]
	end
end