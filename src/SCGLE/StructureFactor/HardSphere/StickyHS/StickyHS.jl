@doc """
`IS_HS_Sticky(τ::Float64, ϕ::Float64, k::Float64)`
Computes the inverse of the static structure factor of a Stucky Hard Sphere.
# Arguments
- `τ::Float64`: .
- `ϕ::Float64`: Volume fraction.
- `k::Float64`: Wave vector.
# References
[1] 

Contributed by O. Joquín'Jaime
"""
function IS_HS_Sticky(τ::Float64, ϕ::Float64, k::Float64)
    # Baxter factorization parameters
    c = (1.0 + 0.5 * ϕ) / (1.0 - ϕ)^2
    b = (ϕ / (1.0 - ϕ)) + τ
    a = ϕ / 12.0

    lm1 = (b + sqrt(b^2 - 4.0 * a * c)) / (2.0 * a)
    lm2 = (b - sqrt(b^2 - 4.0 * a * c)) / (2.0 * a)
    
    λ = min(lm1, lm2)
    σ = 1.0
    μ = λ * ϕ * (1.0 - ϕ)
    
    A = 0.5 * (1.0 + 2.0 * ϕ - μ) / (1.0 - ϕ)^2
    B = σ * (-3.0 * ϕ + μ) / (2.0 * (1.0 - ϕ)^2)
    C = -σ^2 * A - σ * B + σ^2 * λ / 12.0

    cosk = cos(k)
    sink = sin(k)

    α = 1.0 - ((12.0 * ϕ) / k^3) * (2.0 * A * (k * cosk - sink) + (B / σ) * k * (cosk - 1.0) + (λ * k^2 / 12.0) * sink)
    β = (12.0 * ϕ / k^3) * (2.0 * A * (k * sink + cosk - 1.0 - k^2 / 2.0) + (B / σ) * k * (sink - k) + (λ * k^2 / 12.0) * (1.0 - cosk))

    is = α^2 + β^2
    return is
end

@doc """
`S_HS_Sticky(τ::Float64, ϕ::Float64, k::Float64)`
Computes the static structure factor of a Stucky Hard Sphere.
# Arguments
- `τ::Float64`: .
- `ϕ::Float64`: Volume fraction.
- `k::Float64`: Wave vector.
# References
[1] 

Contributed by O. Joquín'Jaime
"""
S_HS_Sticky(τ::Float64, ϕ::Float64, k::Float64) = 1/IS_HS_Sticky(τ, ϕ, k)
