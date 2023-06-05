"""`longtimes(A :: Vector{Float64}, t :: Vector{Float64}, t_long :: Vector{Float64})`
Extends the values of an array following the tendecy of the last points.
# Arguments
- `A :: Vector{Float64}`: Quantity array to extend.
- `t :: Vector{Float64}`: time domain of the Quantity.
- `t_long :: Vector{Float64}`: Extension domain.
"""
function longtimes(A :: Vector{Float64}, t :: Vector{Float64}, t_long :: Vector{Float64})
	# Computing the exponent of growth
	Aⁿ⁺¹ = A[end]
	Aⁿ = A[end-1]
	tⁿ⁺¹ = t[end]
	tⁿ = t[end-1]
	μ = (log(Aⁿ⁺¹)-log(Aⁿ))/(log(tⁿ⁺¹)-log(tⁿ))
	# Computing a second exponent of growth
	Bⁿ⁺¹ = A[end-1]
	Bⁿ = A[end-2]
	τⁿ⁺¹ = t[end-1]
	τⁿ = t[end-2]
	ν = (log(Bⁿ⁺¹)-log(Bⁿ))/(log(τⁿ⁺¹)-log(τⁿ))
	# Selecting the bigger one.
	μ = maximum([μ, ν])
	A_long = (Aⁿ/(tⁿ^μ))*(t_long.^μ)
	# Updating the array
	for i in 1:length(t_long)
		append!(A, A_long[i])
	end
	return A
end
