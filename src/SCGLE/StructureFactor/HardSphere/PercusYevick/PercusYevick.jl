# Percus-Yevick solution's coeficients
α₁(ϕ :: Float64) =  -(1.0 + 2.0*ϕ)^2/(1.0 - ϕ)^4
α₂(ϕ :: Float64) = 6.0*ϕ*(1.0 + 0.5*ϕ)^2/(1.0 - ϕ)^4
α₃(ϕ :: Float64) = 0.5*ϕ*α₁(ϕ)

# Some usefull integrals
# ∫(0,1) dx x² j₀(kx)
I₁(k :: Float64) = k == 0.0 ? 1/3 : (sin(k) - k*cos(k))/k^3
# ∫(0,1) dx x² xj₀(kx)
I₂(k :: Float64) = k == 0.0 ? 1/4 : (-(k^2 - 2.0)*cos(k) + 2.0*k*sin(k) - 2.0)/k^4
# ∫(0,1) dx x² x³j₀(kx)
I₃(k :: Float64) = k == 0.0 ? 1/6 : (4.0*k*(k^2 - 6.0)*sin(k) - (k^4 - 12.0*k^2 + 24.0)*cos(k) + 24.0)/k^6

# FT of Wertheim's direct correlation functions
@doc """
`C_HS_PY(ϕ:: Float64, k :: Float64)` computes the direct correlation function from the Ornstein-Zernike Equation using the Percus-Yevick closure[1].
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
C_HS_PY(ϕ :: Float64, k :: Float64) = α₁(ϕ)*I₁(k) + α₂(ϕ)*I₂(k) + α₃(ϕ)*I₃(k)

@doc """`IS_HS_PY(ϕ:: Float64, k :: Float64)` 
Computes the inverse of the static structure factor using the direct correlation function from the Ornstein-Zernike Equation using the Percus-Yevick closure[1].
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
IS_HS_PY(ϕ :: Float64, k :: Float64) = 1.0 - 24*ϕ*C_HS_PY(ϕ, k)

@doc """`S_HS_PY(ϕ :: Float64, k :: Float64)`
Returns the Static Structure Factor under the Percus-Yevick Closure[1].
```math
S(k) = \\frac{1}{1 - 24\\phi c(k)}
```
# Arguments
- `ϕ :: Float64`: the volume fraction.
- `k :: Float64`: wave vector.
# References
[1] J. P. Hansen and I. McDonald. Theory of Simple Liquids. Academic, London, 1990.
"""
S_HS_PY(ϕ :: Float64, k :: Float64) = 1/IS_HS_PY(ϕ, k)
