"""`adapt(k :: Vector{Float64}, Sf :: Function)`
Updates the wave vector grid at a given waiting time.
       3 +-----------------------------------------------------------------+
         |                +               +            C   +               |
         |                                             a                   |
     2.5 |-+                                          aa                 +-|
         |                                            a                    |
         |                                            aa                   |
       2 |-+                                          a a                +-|
         |                             aAa            a a                  |
         |                            aa aa           a a                  |
     1.5 |-+                        aaa   a           a a                +-|
S(k)     |                        aaa      a          a a                  |
         |                   aaaaaa        a          a a                  |
         |aaaa aaa aaaaaaaaaaa              a         a a  aaa a           |
       1 |-+                                a        a  aaaa aaaaaaaaaa  +-|
         |                                  aa       a   aa                |
         |                                   aa      a                     |
     0.5 |-+                                  aa    aa                   +-|
         |                                     aaBaaa                      |
         |                +               +                +               |
       0 +-----------------------------------------------------------------+
        0.01             0.1              1                10             100
A first maximum
B first minimum
C second maximum
- Between minimum(k) and A a Chebyshev grid
- Between A and C a second Chebyshev gid
# Arguments
- `k :: Vector{Float64}`: Initial wave vector grid.
- `Sf :: Function`: function to evaluate the static structure factor.
"""
function adapt(k :: Vector{Float64}, Sf :: Function)
	S = Sf.(k)
	s = S[k.<5]
	ss = S[k.>5]
	B = findfirst(item-> item == minimum(s), s)
	C = length(s) + findfirst(item-> item == maximum(ss), ss)
	if B == 1 || B == 2
		grid1, dump = ChebyshevGrids(minimum(k), k[C], 200)
		grid2 = collect(k[C]:0.1:maximum(k))
		return unique(vcat(grid1, grid2))
	else
		A = findfirst(item-> item == maximum(s), s)
		if A == 1 || A == 2
			grid1, dump = ChebyshevGrids(minimum(k), k[C], 200)
			grid2 = collect(k[C]:0.1:maximum(k))
			return unique(vcat(grid1, grid2))
		else
			grid1, dump = ChebyshevGrids(minimum(k), k[A], 200)
			grid2, dump = ChebyshevGrids(k[A], k[C],200)
			grid3 = collect(k[C]:0.1:maximum(k))
			return unique(vcat(grid1, grid2, grid3))
		end
	end
end
