"""`blip(T :: Real; ν = 6)`
Returns the effective diameter of an equivalent soft sphere[1] using the fomula
```math
\\lambda^3(T, \\nu) = 1 - 3 \\int_0^1 dx x^2exp(-(1/T)(1/x^(2\\nu)-2/x^\\nu + 1))
```
# Arguments
- `T :: Real`: Temperature.
# Keywords
- `ν = 6`: Parameter to module the softness of the potential.
# References
[1] Luis Enrique Sánchez-Díaz, Pedro Ramírez-González, and Magdaleno Medina-Noyola Phys. Rev. E 87, 052306 – Published 22 May 2013
"""
function blip(T :: Real; ν = 6)
	if T == 0 return 1.0
	else
		dx = 0.001
		x = collect(dx:dx:1)
		integrand = x.*x.*exp.(-(1/T)*(x.^(-2*ν)-2*x.^(-ν).+1))
		λ³ = -3*dx*sum(integrand) + 1
		return λ³^(1/3), λ³
	end
end
