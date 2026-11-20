module PowerElectronics

using DyadData: DyadTimeseries

"""
    channel_resistance(R₀, α₁, α₂, ΔT)

Pristine on-state channel resistance R₀·(1 + α₁·ΔT + α₂·ΔT²), ΔT = T_j − T_ref in K.
The single definition of the channel law: `ConductionMOSFET` and the Stage A components in
`dyad/stageA.dyad` all call this function, so the law cannot differ between them.
Aging enters elsewhere (ΔR in `ConductionMOSFET`, `a` on the thermal side): this is the
pristine law only. A plain polynomial traces symbolically, so Dyad relations can call it.
"""
channel_resistance(R₀, α₁, α₂, ΔT) = R₀ * (1 + α₁ * ΔT + α₂ * ΔT^2)

"""
    stageA_tstops(test)

Solver `tstops` for every burst-resolved solve of device `test`: the time column of
`assets/stageA/test<N>_heat.csv`, which holds every burst edge plus t = 0 and the end pad of
the first segment. The heat input is `ConstantInterpolation` of that series, discontinuous at
exactly these points, and bursts can be 35 ms apart, so a solver that does not land on the
edges steps over whole bursts. Needs the data files (see the README).
"""
stageA_tstops(test::Integer) =
    Float64.(DyadTimeseries("dyad://PowerElectronics/stageA/test$(test)_heat.csv";
                            independent_var = "time", dependent_vars = ["P"])["time"])

# Observation model: seeded measurement noise on saved samples
include("observation.jl")

include("../generated/module.jl")

end # module PowerElectronics
