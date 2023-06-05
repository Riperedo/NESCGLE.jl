"""`ChebyshevGrids(a :: Float64, b :: Float64, N :: Integer)`
Returns a Chevishev nodes and weights[1]
# Arguments
- `a :: Float64`: Initial node.
- `b :: Float64`: Final node.
- `N :: Integer`: Number of nodes.
# References
[1]  C.W. Clenshaw, A.R. Curtis, A method for numerical integration on an automatic computer, Numer. Math. 2 (1960) 197–205.
"""
function ChebyshevGrids(a :: Float64, b :: Float64, N :: Integer)
	@assert a < b "The first argument must be less than the second."
	i = collect(1:N)
	x = map(i -> 0.5*(a+b) + 0.5*(b-a)*cospi((2*i-1)/(2*N)), i)
	y = [-1.0]
	append!(y, map(i -> cospi((2*i-1)/(2*N)), i))
	append!(y, 1.0)
	w = map(y -> π*0.5*(b-a)*sqrt(1.0 - y*y)/N, y)
	append!(x, a)
	append!(x, b)
	sort!(x)
	return x, w
end

@doc """UniformGrid(a :: Float64, b :: Float64, step :: Float64) = collect(a:step:b)
Returns a uniform grid. Is exactly the same as the fucntion `collect`.
# Arguments
# Arguments
- `a :: Float64`: Initial point.
- `b :: Float64`: Final point.
- `step :: Float64`: width step.
"""
UniformGrid(a :: Float64, b :: Float64, step :: Float64) = collect(a:step:b)

"""`ExponentialGrid(a :: Int64, b :: Int64, x₀ :: Float64)`
Return an exponential grid of the form
```math
X = {x vert x = x_010^{n/a}text{ for } n in [0,b]}
```
the sugested parameters are a = 5, b = 60.
# Arguments
- `a :: Int64`: Number por points between decade.
- `b :: Int64`: Number of elements in the array.
- `x₀ :: Float64`: Initial value.
"""
function ExponentialGrid(a :: Int64, b :: Int64, x₀ :: Float64)
	@assert x₀ > 0.0 "x₀ must be positive"
	return x₀*[10^(n/a) for n in 0:b]
end
