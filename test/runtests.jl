using PowerElectronics
using Test

# The NASA-derived input files in assets/stageA are not distributed with the package (see the
# README). Tests that read them run only when the files are present.
const HAS_DATA = isdir(joinpath(pkgdir(PowerElectronics), "assets", "stageA"))

include("../generated/tests.jl")
include("components.jl")      # component-level units: fail here → one component, usually one term
include("physics.jl")         # system-level physics of the rigs (KVL, energy balance, steady state)
include("thermal.jl")         # CaseToAmbient law, energy, τ growth, reversed gradient, radiation
include("controller.jl")      # hysteresis controller units + closed-loop smoke (test 30)
include("observation.jl")     # observation blocks CurrentOffset, VdsStepArtefact, ThermocoupleLag; noise helper
include("aging.jl")

if HAS_DATA
    include("stageA_data.jl")          # data contracts of the Stage A input files (assets/stageA)
    include("stageA_model.jl")         # Stage A model: energy balance, junction steady state
    include("observation_harness.jl")  # ObservationHarness
    include("stageA_lag.jl")           # StageAModelLag (StageAModel + ThermocoupleLag)
else
    @info "assets/stageA not found: skipping the tests that need the NASA-derived data files"
end
