# CaseTemperatureController unit checks.
#
# The gate g and latch L are event-updated PARAMETERS; their history is not tracked by the
# solution (sol(t) evaluates parameter-dependent observables at the final parameter state), so
# the assertions here are staged FINAL-state checks on short horizons plus a closed-loop smoke
# on the saved T_c trajectory — a true state, which the solution does track.

using ModelingToolkit, DyadInterface, OrdinaryDiffEqDefault
using Statistics: mean, quantile
const getp = ModelingToolkit.getp   # re-exported SII getter; SII itself is not a test dep

"Final (g, L) of the controller after `stop` seconds at a constant measured temperature."
function ctl_final(Tc_const, stop)
    c = PowerElectronics.CaseTemperatureController(name = :c)   # band 483.15/485.15, shutdown 623.15, T_sample 2.5
    vars = ModelingToolkit.@variables d(ModelingToolkit.t_nounits)
    w = ModelingToolkit.System([c.Tc ~ Tc_const, vars[1] ~ c.enable],
                               ModelingToolkit.t_nounits; systems = [c], name = :w)
    sys = DyadInterface.simplify_model(w)
    prob = ODEProblem(sys, Dict{Any, Any}(), (0.0, stop))
    sol = solve(prob)
    (g = getp(sys, c.g)(sol), L = getp(sys, c.L)(sol), retcode = sol.retcode)
end

@testset "CaseTemperatureController" begin
    r = ctl_final(480.0, 6.0)      # below the band: stays on
    @test string(r.retcode) == "Success"
    @test r.g == 1.0 && r.L == 1.0
    r = ctl_final(500.0, 6.0)      # above hi, below shutdown: off, not latched
    @test r.g == 0.0 && r.L == 1.0
    r = ctl_final(650.0, 6.0)      # above shutdown: latched off
    @test r.L == 0.0 && r.g == 0.0
    # latch permanence: the ControllerLatch sine ends far below the band at t = 400 (270 K);
    # a broken latch would re-arm there
    rl = PowerElectronics.ControllerLatch()
    ml = DyadInterface.symbolic_container(rl)
    @test string(rl.sol.retcode) == "Success"
    @test getp(ml, ml.ctl.L)(rl.sol) == 0.0
    @test getp(ml, ml.ctl.g)(rl.sol) == 0.0
    # re-arm: the ControllerCycle sine ends below the band — the gate must be back on
    rc = PowerElectronics.ControllerCycle()
    mc = DyadInterface.symbolic_container(rc)
    @test getp(mc, mc.ctl.g)(rc.sol) == 1.0
    @test getp(mc, mc.ctl.L)(rc.sol) == 1.0
end

@testset "ControlledAveragedRig closed loop, test 30" begin
    r = PowerElectronics.ControlledRig30()
    m = DyadInterface.symbolic_container(r)
    @test string(r.sol.retcode) == "Success"
    tcm = [r.sol(x; idxs = m.T_c_meas) for x in 400.0:1.0:1741.0] .- 273.15
    # regulates to the 209–210 °C band (measured hold mean 208.4 °C) and never trips the
    # 220 °C shutdown — the latch would freeze T_c into a monotone cool-down
    @test 205.0 <= mean(tcm) <= 212.0
    @test maximum(tcm) < 219.0
    # it cycles: the hold is not a flat line (measured p2–p98 ripple 12.8 °C; model ~15)
    rip = quantile(tcm, 0.98) - quantile(tcm, 0.02)
    @test 8.0 <= rip <= 20.0
end
