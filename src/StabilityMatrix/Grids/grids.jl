"""
    ChebyshevGrids(a :: Float64, b :: Float64, N :: Integer)
    
Returns a Chebyshev nodes and weights[1]
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

@doc"""
    UniformGrid(a :: Float64, b :: Float64, step :: Float64) = collect(a:step:b)
Returns a uniform grid. Is exactly the same as the fucntion `collect`.
# Arguments
- `a :: Float64`: Initial point.
- `b :: Float64`: Final point.
- `step :: Float64`: width step.
"""
UniformGrid(a :: Float64, b :: Float64, step :: Float64) = collect(a:step:b)

"""
    ExponentialGrid(a :: Int64, b :: Int64, x₀ :: Float64)
Return an exponential grid of the form
`X = {x | x = x10^(n/a) for x in [0,b]}`
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

#############
#	Grids	#
#############

"""
    Cheb_plus_U(x₀ :: Float64, x₁ :: Float64, x₂ :: Float64, N₁ :: Int64, N₂ :: Int64)
Auxiliar function that returns a composed grid by a Chebishev grid and a uniform grid.
# Arguments
- `x₀ :: Float64`: Initial point.
- `x₁ :: Float64`: Middle point.
- `x₂ :: Float64`: Final porint.
- `N₁ :: Int64`: Number of point in the Chebishev grid.
- `N₂ :: Int64`: number of point in the uniform grid.
"""
function Cheb_plus_U(x₀ :: Float64, x₁ :: Float64, x₂ :: Float64, N₁ :: Int64, N₂ :: Int64)
	# first grid
	Chev, w = ChebyshevGrids(x₀, x₁, N₁)
	x = vcat(Chev, collect(x₁:((x₂-x₁)/N₂):x₂))
	return unique!(x)
end

"""
    grid_plus_grid(x :: Array{Float64}, y :: Array{Float64})
Auxiliar fucntion to concat two grids.
# Arguments
- `x :: Array{Float64}`: First grid.
- `y :: Array{Float64}`: Second grid.
"""
function grid_plus_grid(x :: Array{Float64}, y :: Array{Float64}) 
	@assert x[end] <= y[1] "The first grid must have a domain lower than the second one. 
	The last element of the first grid is "*string(x[end])*", and the first element for the second array is "*string(y[1])*"."
	return unique!(vcat(x, y))
end
