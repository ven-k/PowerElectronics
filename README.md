# PowerElectronics

A Dyad component library for an electro-thermal power MOSFET (IRF520N) with an aging state and an explicit
model of its sensors. It compiles to ModelingToolkit. The component values were identified from NASA's MOSFET
thermal-overstress run-to-failure data ([Data](#data)).

The library is under the BSD 3-Clause licence ([LICENSE](LICENSE)). It runs on Dyad, which is free of charge for
non-commercial use only, and on ModelingToolkit, whose tearing dependencies are AGPL-3.0.
[Licence terms](#licence-terms) lists them.

What the library does not contain: a validated law for how the aging state grows. The laws we fitted did not
pass their validation, so `AgingState` ships without one and you supply your own.

## Components

Every component's docstring states what it models and where each default value comes from.

| Component | File | What it models |
|---|---|---|
| `ConductionMOSFET` | `dyad/mosfet.dyad` | On-state channel resistance R₀·(1 + α₁·ΔT_j + α₂·ΔT_j²) + ΔR, off-state leakage, and the conduction loss as heat into the junction. Switching edges are not modelled. |
| `CauerNetwork` | `dyad/thermal.dyad` | Junction-to-case thermal ladder of `N` RC stages. The aging state `a` raises the resistance of the first stage only: R_th[1]·(1 + a). |
| `DutyAveragedHeat` | `dyad/thermal.dyad` | Cycle-averaged heat of a PWM-driven device: the network receives duty × the on-state loss. |
| `CaseToAmbient` | `dyad/thermal.dyad` | Case-to-ambient path of an unsinked TO-220: free convection (heat flow ∝ ΔT^1.25) plus radiation. |
| `JunctionStage` | `dyad/stageA.dyad` | One junction RC stage with scalar parameters, R_jc·(1 + a) and C_die. |
| `StageAModel` | `dyad/stageA.dyad` | One device as fitted in Stage A: measured burst heat into a junction stage, a case node, a fixture side node and `CaseToAmbient`. |
| `StageAJunction`, `StageAJunction2` | `dyad/stageA.dyad` | Junction (one or two stages) and channel law with the case pinned to the measured case temperature. |
| `AgingState` | `dyad/aging.dyad` | The aging state a(t) with its start value and a junction-temperature input. Partial: you add the growth law. |
| `PWMSource` | `dyad/sources.dyad` | Rectangular gate command (1 kHz, 40 % duty by default). |
| `NASARig` | `dyad/rig.dyad` | The NASA test circuit with a switching gate, for millisecond windows. |
| `AveragedRig` | `dyad/rig.dyad` | The same circuit with cycle-averaged heat, for seconds to hours. |
| `CaseTemperatureController` | `dyad/controller.dyad` | The rig's hysteresis controller on the case temperature, with a shutdown latch. |
| `ControlledAveragedRig`, `Rig30`, `Rig36` | `dyad/controller.dyad` | The averaged rig in closed loop with the controller; `Rig30` and `Rig36` carry the fitted values and logged schedules of tests 30 and 36. |
| `CurrentOffset` | `dyad/observation.dyad` | Zero offset of the drain-current channel. |
| `VdsStepArtefact` | `dyad/observation.dyad` | Step-response artefact of the drain-source voltage channel. |
| `ThermocoupleLag` | `dyad/observation.dyad` | First-order lag of the thermocouple reading of the case temperature. |
| `StageAModelLag` | `dyad/observation.dyad` | `StageAModel` with `ThermocoupleLag` on its case-temperature output. |
| `ObservationHarness` | `dyad/observation.dyad` | Junction and channel observed through `CurrentOffset` and `VdsStepArtefact`. |

Components whose names end in `Test` are fixtures of the test suite.

## Aging

`CauerNetwork`, `JunctionStage` and `ConductionMOSFET` take their aging term as a parameter (`a`, `a`, `ΔR`), which
stays constant during a simulation. Built with `aging_input = true`, each takes it from a signal input instead
(`a_in`, `a_in`, `ΔR_in`), so the term can change while the simulation runs.

`AgingState` holds the state. It has no growth law. Extend it and add one relation:

```
component MyAging
  extends AgingState
  parameter k::Real = 1e-6
relations
  der(a) = k * (a + 0.004)
end
```

Connect a temperature sensor on the junction to `T_j`, and `a` to the `a_in` of the thermal component. For the
electrical term, connect `a` through a gain to `ΔR_in`. `AgingClosedFormTest` in `dyad/test_harnesses.dyad` is a
complete example with a test law; that law is not validated.

## Current zero offset

`offset_I` of `CurrentOffset` and `ObservationHarness` defaults to 0 A, which leaves the reading unchanged. Set it
to the test's measured zero, the median off-state drain-current reading: 0.2644 A for test 30, 0.2480 A for test 36.

Do not fit it together with the channel. With the offset free, the fit returns 2.50 A on test 30 and the
temperature coefficient α₁ is recovered only within 41.5 %; with the offset fixed at the measured zero, within
14.4 %.

## Installation

1. Install Dyad 3.5 by following the Dyad documentation: <https://help.juliahub.com/dyad/dev/>. Dyad is free for
   non-commercial use; see [Licence terms](#licence-terms).
2. Clone this repository and compile the Dyad sources. The repository holds only the `.dyad` sources; the
   compiler writes the Julia code into `generated/`.

   ```
   git clone https://github.com/ven-k/PowerElectronics
   cd PowerElectronics
   dyad compile .
   ```
3. Add it to an environment of your own with `Pkg.develop(path = "path/to/PowerElectronics")`, or work in the
   project itself with `julia --project=.` and `using Pkg; Pkg.instantiate()`.

Run `dyad compile .` again after editing any `.dyad` file. The compiler calls `julia` from `PATH`, which must be
Dyad's Julia build. The package is not registered in the General registry.

## Data

The input files in `assets/stageA/` (per-device burst heat and case-temperature series, derived from the NASA
data) are not distributed with this repository. The components that read them are `StageAModel`,
`StageAJunction`, `StageAJunction2`, `StageAModelLag` and `ObservationHarness`, and the analysis `StageACase30`;
all other components run without them. To obtain the files, please reach out by opening an issue on this
repository.

The raw data are NASA's: J. R. Celaya, A. Saxena, S. Saha and K. Goebel, "MOSFET Thermal Overstress Aging Data
Set," NASA Prognostics Data Repository, NASA Ames Research Center, Moffett Field, CA.
<https://www.nasa.gov/intelligent-systems-division/discovery-and-systems-health/pcoe/pcoe-data-set-repository/>

## Examples

One example per component. Each runs in under 15 s on first use. Run them in the compiled project
(`julia --project=.`).

```julia
using PowerElectronics, DyadInterface, ModelingToolkit, OrdinaryDiffEqDefault
using DyadInterface: simplify_model

# ConductionMOSFET: on-state resistance 100 K above T_ref with ΔR = 0.05 Ω
res = MOSFETConduction(); m = symbolic_container(res)
res.sol[m.mosfet.R_on][end]                      # 0.4766 Ω

# CauerNetwork: junction node after a 2 W heat step into a 3.1 K/W ladder, case at 300 K
res = CauerStep(); m = symbolic_container(res)
res.sol[m.jc.T[1]][end]                          # 306.2 K

# DutyAveragedHeat: a heat capacity charged with duty × the device heat
res = DutyAveraging(); m = symbolic_container(res)
res.sol[m.store.T][end]                          # 310.0 K

# CaseToAmbient: free cool-down from 250 °C over 1200 s
res = CaseCooldown(); m = symbolic_container(res)
res.sol[m.T_c][end]                              # 306.06 K

# JunctionStage: 2 W into 3.1 K/W, case at 400 K
res = JunctionStep(); m = symbolic_container(res)
res.sol[m.jc.T_j][end]                           # 406.2 K

# AgingState: a test law grows the state from 0.1 over 100 s
res = AgingClosedForm(); m = symbolic_container(res)
res.sol[m.law.a][end]                            # 0.7685

# PWMSource: low and high level of the gate command
res = PWMWindow(); m = symbolic_container(res)
extrema(res.sol[m.pwm.y])                        # (0.5, 2.5)

# NASARig: drain current in the middle of the first on-pulse
res = NASARigSwitching(); m = symbolic_container(res)
res.sol(0.2e-3, idxs = m.mosfet.i)               # 16.88 A

# AveragedRig: case temperature after 600 s in open loop at 40 % duty
res = AveragedRigWarmup(); m = symbolic_container(res)
res.sol[m.T_c][end]                              # 692.1 K

# CaseTemperatureController: the shutdown latch after the case crossed shutdownTemp
res = PowerElectronics.ControllerLatch(); m = symbolic_container(res)
ModelingToolkit.getp(m, m.ctl.L)(res.sol)        # 0.0 (latched off)

# ControlledAveragedRig, Rig30 and Rig36: thermocouple reading at the end of the closed-loop window
res = PowerElectronics.ControlledRig30(); m = symbolic_container(res)
res.sol[m.T_c_meas][end]                         # 477.09 K
res = PowerElectronics.ControlledRig36(); m = symbolic_container(res)
res.sol[m.T_c_meas][end]                         # 518.77 K

# CurrentOffset: reading minus input
res = CurrentOffsetRamp(); m = symbolic_container(res)
res.sol[m.meas.y][end] - res.sol[m.meas.u][end]  # 0.3 A

# VdsStepArtefact: largest artefact error after a voltage step
res = VdsArtefactStep(); m = symbolic_container(res)
maximum(abs.(res.sol[m.art.e]))                  # 0.155 V

# ThermocoupleLag: reading 5 s after a step of the case temperature
res = LagStep(); m = symbolic_container(res)
res.sol[m.lag.y][end]                            # 493.15 K
```

These need the data files in `assets/stageA/` ([Data](#data)):

```julia
# StageAModel: test 30's warm-up window driven by its measured burst heat
res = StageATest30(); m = symbolic_container(res)
res.sol[m.T_c][end]                              # 484.02 K

# StageAJunction: channel resistance with the case pinned to test 30's measured temperature
edges = filter(<=(400.0), PowerElectronics.stageA_tstops(30))
sys = simplify_model(PowerElectronics.StageAJunction(; name = :junction))
sol = solve(ODEProblem(sys, [], (0.0, 400.0)); tstops = edges)
sol[sys.R_model][end]                            # 0.706 Ω

# StageAModelLag: case temperature minus thermocouple reading at 400 s
sys = simplify_model(PowerElectronics.StageAModelLag(; name = :lagged))
sol = solve(ODEProblem(sys, [], (0.0, 400.0)); tstops = edges)
sol[sys.T_c][end] - sol[sys.T_c_meas][end]       # 0.196 K

# ObservationHarness: measured minus true drain current with test 30's measured zero offset
sys = simplify_model(PowerElectronics.ObservationHarness(; name = :harness))
sol = solve(ODEProblem(sys, [sys.offset_I => 0.2644], (0.0, 400.0)); tstops = edges)
sol[sys.I_d_meas][end] - sol[sys.I_d][end]       # 0.2644 A
```

## Tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Compile first (`dyad compile .`). The tests that read the data files run only when `assets/stageA/` is present.
There is no hosted continuous integration, because the Dyad packages need a Dyad login to install.

## Licence terms

Our code is under the BSD 3-Clause licence ([LICENSE](LICENSE)). The repository does not contain Dyad's compiler
output or any Dyad package code. The packages the library depends on have their own terms:

| Package (version tested) | Licence |
|---|---|
| DyadInterface 7.5.1, DyadData 2.3.0, DyadModelOptimizer 15.5.1 | JuliaHub's Dyad terms: free of charge for non-commercial use; commercial use needs a licence from JuliaHub. Under those terms, use by a government body counts as commercial. |
| ElectricalComponents, ThermalComponents, RotationalComponents, TranslationalComponents 3.0.0; BlockComponents 5.0.0 | BSD 3-Clause |
| ModelingToolkit 11.45.3 | MIT |
| ModelingToolkitTearing 1.20.7, StateSelection 1.11.1 (dependencies of ModelingToolkit); DyadCompilerPasses 0.3.7 | AGPL-3.0 |
