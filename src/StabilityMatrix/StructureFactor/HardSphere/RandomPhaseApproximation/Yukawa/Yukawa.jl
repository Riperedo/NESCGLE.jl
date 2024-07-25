"""`βU_Yukawa(A :: Float64, z :: Float64, k :: Float64)`
Auxiliar function to compute the Fourier Transform of Yukawa potential[1].
```math
\\beta u(k) = -4\\pi(A\\sigma^3)\\frac{k\\sigma\\cos(k\\sigma)+z\\sin(k\\sigma)}{k\\sigma[(k\\sigma)^2+z^2]}
```
# Arguments
- `A :: Float64`: Aplitude potential (effective charge).
- `z :: Float64`: Inverse screening length.
- `k :: Float64`: Wave vector.
# References
[1] Yukawa, H. (1935). "On the interaction of elementary particles". Proc. Phys.-Math. Soc. Jpn. 17: 48.
"""
function βU_Yukawa(A :: Float64, z :: Float64, k :: Float64)
	if k < 0.075
		return 4*π*A*((1/z)+1-(z*z+9*z+6)*k*k/(6*z*z))
	else
		return 4*π*A*(k*cos(k) + z*sin(k))/(k*(k^2+z^2))
	end
end
