# Observation model, the noise part: T_c,meas = T_c + noise.
#
# The Dyad blocks of dyad/observation.dyad are deterministic. Noise is added here, to saved samples, for two
# reasons: BlockComponents 5.0.0 has no noise source, and a deterministic ODE keeps any gradient path through
# the model free of random terms.

using Random: Xoshiro, randn

"""
    add_measurement_noise(x, σ; seed) -> Vector{Float64}

Return `x .+ σ .* ε`, with ε independent standard normal draws from `Random.Xoshiro(seed)`: the
observation noise added to saved samples `x` (e.g. a `ThermocoupleLag` reading in K, with σ = σ_Tc,
the standard deviation of the case-temperature noise).

- `seed` is required, so every noisy series is reproducible; the same `seed` and length give the same
  draws within one Julia version (the stream of `Xoshiro` may change between Julia versions).
- `x` is not modified. `σ = 0` returns a copy of `x` as `Float64`.
- `σ` must be finite and non-negative.
"""
function add_measurement_noise(x::AbstractVector{<:Real}, σ::Real; seed::Integer)
    (isfinite(σ) && σ >= 0) || throw(ArgumentError("σ must be finite and ≥ 0, got $σ"))
    ε = randn(Xoshiro(seed), length(x))
    return Float64.(x) .+ Float64(σ) .* ε
end
