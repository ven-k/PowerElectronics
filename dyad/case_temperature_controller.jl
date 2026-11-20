# Implementation of the external component `CaseTemperatureController` (dyad/controller.dyad).
#
# Hysteresis with memory needs discrete state, which Dyad has no construct for. Here the gate g
# and the shutdown latch L are parameters held constant between samples and updated by a
# periodic ModelingToolkit discrete event every `T_sample` seconds, as the rig's digital
# controller does.
#
# The setpoint schedule is piecewise-constant in time (the warm-up staircase advances
# on a timed dwell). It enters as plain vectors at construction — structural, like a dataset —
# through two opaque callable parameters bound to `StairSchedule` lookups (the same callable-
# parameter pattern the stdlib Interpolation block and the KinematicPTPBoundedJerk example use).
#
# At each sample, Tc ≥ hi commands off, Tc ≤ lo commands on, and Tc ≥ shutdownTemp zeroes the
# latch permanently (nothing ever sets L back to 1), so no up/down logic is needed.

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, Pre

"Piecewise-constant schedule: value of the step whose start time is the last one ≤ x."
struct StairSchedule
    ts::Vector{Float64}
    vs::Vector{Float64}
end
(s::StairSchedule)(x) = s.vs[clamp(searchsortedlast(s.ts, x), 1, length(s.vs))]

function CaseTemperatureController(; name, N = 1, t_sched = [0.0], lo_sched = [483.15],
        hi_sched = [485.15], T_sample = 2.5, shutdownTemp = 623.15, duty_on = 0.4)
    length(t_sched) == length(lo_sched) == length(hi_sched) ||
        error("CaseTemperatureController: schedule vectors must share one length")
    vars = @variables begin
        (Tc(t)::Real), [input = true, description = "Measured case temperature (K)"]
        (duty_out(t)::Real), [output = true, description = "Duty command (0 while gated off)"]
        (enable(t)::Real), [output = true, description = "Gate enable (0/1), latched 0 after shutdown"]
        (lo_now(t)::Real), [description = "Active lower setpoint (schedule lookup)"]
        (hi_now(t)::Real), [description = "Active upper setpoint (schedule lookup)"]
    end
    pars = @parameters begin
        shutdownTemp = shutdownTemp, [description = "Latch-off temperature (K)"]
        duty_on = duty_on, [description = "Commanded duty while enabled"]
        # Discrete controller states as PARAMETERS with NUMERIC defaults: a symbolic default
        # (g = gate0) becomes an immutable BINDING that affect writes cannot touch, and a
        # D(x) ~ 0 state is constant-folded away. Plain parameters are affect-writable.
        g = 1.0, [description = "Hysteresis gate state (updated by the sampler)"]
        L = 1.0, [description = "Shutdown latch (updated by the sampler)"]
    end
    lo_f = only(@parameters (lo_f::StairSchedule)(..))
    hi_f = only(@parameters (hi_f::StairSchedule)(..))
    bindings = Dict(
        lo_f => StairSchedule(Float64.(collect(t_sched)), Float64.(collect(lo_sched))),
        hi_f => StairSchedule(Float64.(collect(t_sched)), Float64.(collect(hi_sched))))
    Tc, duty_out, enable, lo_now, hi_now = vars
    eqs = [
        enable ~ pars[3] * pars[4],
        duty_out ~ pars[2] * pars[3] * pars[4],
        lo_now ~ lo_f(t),
        hi_now ~ hi_f(t),
    ]
    # SAMPLED hysteresis: the real rig's controller is a digital loop — the log cadence equals
    # rows/duration (2.48 s test 30 run 1, 1.41 s test 36) and the measured hold period is ~2
    # samples with ~1 sample on-time. A periodic discrete event reads the lagged Tc every
    # T_sample and updates the held states. Continuous crossing detection both over-toggles at
    # the fast early staircase steps (solver chatter to MaxIters) and is the wrong model of a
    # sampled controller. T_sample is structural: an event period must be concrete.
    # The affect is IMPERATIVE (plain Julia on numeric values): in this MTKBase/SymbolicUtils
    # pairing, symbolic equation-affects break twice over — self-references compile into an
    # implicit discrete problem, and Pre(x) terms carry Unknown shape that fails
    # promote_shape on any comparison rebuild. ImperativeAffect reads the lagged Tc and the
    # active band, updates the held states, and returns them — no symbolic codegen at all.
    sampler = ModelingToolkit.ImperativeAffect(
        (m, o, c, i) -> begin
            Lnew = o.Tc >= o.sd ? 0.0 : m.L
            gnew = (o.Tc >= o.hi ? 0.0 : (o.Tc <= o.lo ? 1.0 : m.g)) * Lnew
            (; g = gnew, L = Lnew)
        end;
        observed = (; Tc = Tc, lo = lo_now, hi = hi_now, sd = pars[1]),
        modified = (; g = pars[3], L = pars[4]))
    evs_d = [Float64(T_sample) => sampler]
    System(eqs, t, vars, [pars; lo_f; hi_f]; name, discrete_events = evs_d, bindings)
end
