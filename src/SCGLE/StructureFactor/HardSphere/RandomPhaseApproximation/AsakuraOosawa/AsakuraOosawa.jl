"""`βU_AO_1(ϕ :: Real, nₚ :: Real, ξ :: Real, k :: Real)`
Monodisperse approximation of Asakura Oosawa potential[1]
where the AO model is replaced by a Yukaya potential
```math
b = 3n_p(1+\\xi)/\\xi^2 \\quad A = exp[3n_p(1+\\xi)/(2\\xi)]
```
where nₚ is the density of the polymer and ξ=2R/σ with R is the polymer’s radius of gyration.
# References
[1] J. Bergenholtz, W. C. K. Poon, and M. Fuchs Langmuir 2003 19 (10), 4493-4503 DOI: 10.1021/la0340089
"""
function βU_AO_1(ϕ :: Real, nₚ :: Real, ξ :: Real, k :: Real)
	b = 3*nₚ*(1+ξ)/(ξ^2)
	A = exp(3*nₚ*(1+ξ)/(2*ξ))
	return -4*π*A*(k*cos(k) + b*sin(k))/(k*(k^2+b^2))
end

""" `V_overlap_AO(q :: Real, ξ :: Real)`
Auxiliary function to evaluate the volume overlap in the Asakura Oosawa Model.
Vₒᵥₑᵣₗₐₚ the so-called overlap volume given by
```math
V_\\text{overlap}= (\\pi/6)\\sigma^3(1+\\xi)^3[1 - (3/2)(r/\\sigma)(1+\\xi)^{-1} + (1/2)(r/\\sigma)^3(1+\\xi)^{-3}]
```
# Arguments
- `q :: Float64`: wave vector magnitude.
- `ξ :: Float64`: polimers' radius of gyration.
"""
function V_overlap_AO(q :: Real, ξ :: Real)
	a = 1 + ξ
	aq = a*q
	cosaq = cos(aq)
	sinaq = sin(aq)
	cosq = cos(q)
	sinq = sin(q)
	V = 0
	if q > 0.1
		A = (sinaq-aq*cosaq-sinq+q*cosq)/(q^3)
		B = ((2-aq^2)*cosaq-2*q*(sinq-a*sinaq) + (q^2-2)*cosq)/(a*(q^4))
		C = (4*aq*(aq^2-6)*sinaq - (aq^4-12*(aq^2)+24)*cosaq- 4*q*(q^2-6)*sinq + (q^4-12*(q^2)+24)*cosq)/((a^3)*(q^6))
		V = (A - (3/2)*B + (1/2)*C)
	else
		A = (1/3)*(a^3 - 1) + (1/30)*(1 - a^5)*q^2 + (1/840)*(a^7 - 1)*q^4
		B = (a^4 - 1)/(4*a) - ((a^6 - 1)*q^2)/(36*a) + ((a^8 - 1)*q^4)/(960*a)
		C = (a^6 - 1)/(6*a^3) - ((a^8 - 1)*q^2)/(48*a^3) + ((a^10 - 1)*q^4)/(1200*a^3) 		
		V = (A - (3/2)*B + (1/2)*C)
	end
	return -4*π*(π/6)*(a^3)*V
end

"""`βU_AO_2(ϕ :: Real, nₚ⁽ᴿ⁾ :: Real, ξ :: Real, k :: Real)`
Monodisperse approximation of Asakura Oosawa potential[1] where the AO model is replaced by a Square Well potential
```math
\\beta U(r) = \\begin{cases}
\\infty & r<\\sigma \\\\
-\\beta\\Pi_pV_\\text{overlap} & r\\in[\\sigma, \\sigma(1+\\xi)]\\\\
0 & \\sigma(1+\\xi) < r
\\end{cases}
```
where Πₚ is the osmotic pressure of a polymer solution in equilibrium, given by βΠₚ = nₚ⁽ᴿ⁾ and nₚ⁽ᴿ⁾ must be interpreted as the number of polymer coils in the volume accessible to them.
nₚ = αnₚ⁽ᴿ⁾[2]
where
α = (1-ϕ)exp(-Aγ-Bγ²-Cγ³)
and
γ = ϕ(1-ϕ)⁻¹
A = 3ξ + 3ξ² + ξ³
B = 9ξ²/2 + 3ξ³
C = 3ξ³
nₚ is the density of the polymer and ξ=2R/σ with R is the polymer’s radius of gyration.
# References
[1] J. Bergenholtz, W. C. K. Poon, and M. Fuchs Langmuir 2003 19 (10), 4493-4503 DOI: 10.1021/la0340089
[2] H. N. W. Lekkerkerker et al 1992 EPL 20 559
"""
function βU_AO_2(ϕ :: Real, nₚ⁽ᴿ⁾ :: Real, ξ :: Real, k :: Real)
	#a = 1 + ξ
	#γ = ϕ/(1-ϕ)
	#A = 3*ξ + 3*(ξ^2) + ξ^3
	#B = 9*(ξ^2)/2 + 3*(ξ^3)
	#C = 3*(ξ^3)
	#α = (1-ϕ)*exp(-A*γ-B*(γ^2)-C*(γ^3))
	#Πₚ = nₚ/α
	Πₚ = nₚ⁽ᴿ⁾
	return -Πₚ*V_overlap_AO(k, ξ)
end
