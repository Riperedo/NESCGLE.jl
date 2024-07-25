@doc """ϕ_VW(ϕ :: Float64)
Comutes the density correction for hard spheres using the Verlet-Weiss approximation[1].
# References
[1] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
"""
ϕ_VW(ϕ :: Float64) = ϕ*(1.0 - (ϕ / 16.0))

@doc """k_VW(ϕ :: Float64, k :: Float64)
Comutes the wave vector correction for hard spheres using the Verlet-Weiss approximation[1].
# References
[1] Loup Verlet and Jean-Jacques Weis. Phys. Rev. A 5, 939 – Published 1 February 1972
"""
k_VW(ϕ :: Float64, k :: Float64) = k*((ϕ_VW(ϕ)/ϕ)^(1.0/3.0))
