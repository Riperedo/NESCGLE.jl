using ModeCouplingTheory
using StaticArrays, LinearAlgebra

import ModeCouplingTheory.evaluate_kernel
import ModeCouplingTheory.evaluate_kernel!
import ModeCouplingTheory.MemoryKernel

#############
#   MSD     #
#############
struct MSD_SCGLEKernel{TDICT,T} <: MemoryKernel
    tDict::TDICT
    Δζ::T
end

"""
MSDS_CGLEKernel(ϕ, D⁰, k_array, Sₖ, sol)

Constructor of a MSD_SCGLEKernel. It implements the kernel

K(k,t) = D⁰ / (36πϕ) ∫dq q^4 c(q)^2 F(q,t) Fs(q,t)

where the integration runs from 0 to infinity. F and Fs are the coherent
and incoherent intermediate scattering functions, and must be passed in
as solutions of the corresponding equations.

# Arguments:

* `ϕ`: volume fraction
* `D⁰`: Short times diffusion coefficient
* `k_array`: vector of wavenumbers at which the structure factor is known
* `Sₖ`: vector with the elements of the structure factor 
* `sol`: a solution object of an equation with a SCGLEKernel.

# Returns:

an instance `k` of `MSD_SCGLEKernel <: MemoryKernel`, which can be evaluated like:
`k = evaluate_kernel(kernel, F, t)`
"""
function MSD_SCGLEKernel_constructor(sol, k_array)
    Δζ = []
    t = sol.t
    Nt = length(t)
    kc = 1.305*2π
    λ = 1/(1+(k_array[1]/kc)^2)
    K = get_K(sol)
    for i in 1:Nt
        append!(Δζ, K[i][1,1]/λ)
    end
    tDict = Dict(zip(t, eachindex(t)))
    kernel = MSD_SCGLEKernel(tDict, Δζ)
    return kernel
end

function evaluate_kernel(kernel::MSD_SCGLEKernel, MSD, t)
    it = kernel.tDict[t]
    return kernel.Δζ[it]
end

function MSD(k::Vector{Float64}, sol, solver)
    MSD0 = 0.0; dMSD0 = 0.0; α = 0.0; β = 1.0; γ = 0.0; δ = -1.0;
    msdkernel = MSD_SCGLEKernel_constructor(sol, k)
    msdequation = MemoryEquation(α, β, γ, δ, MSD0, dMSD0, msdkernel)
    msdsol = solve(msdequation, solver)
    return msdsol.t, msdsol.F, msdkernel.Δζ
end

@doc"""
ΔG = (kBT/60π²)∫dk k⁴[(1/S)(∂S/∂k)]²[(F/S)]^2
"""
function get_ΔG(sol, k_array::Vector{Float64}, Sₖ::Vector{Float64})
    t = sol.t
    Nk = div(length(k_array),2)
    Δk = k_array[2] - k_array[1]
    k⁴ = k_array[Nk+1:end-1].^4
    ∂S = diff(Sₖ[Nk+1:end])./diff(k_array[Nk+1:end])
    S = Sₖ[Nk+1:end-1]
    F = sol.F
    ΔG = []
    for i in 1:length(t)
        dG = sum(k⁴.*((∂S./S).^2).*((F[i][Nk+1:end-1]./S).^2))
        dG *= Δk/(60*π*π)
        append!(ΔG, dG)
    end
    return ΔG
end

