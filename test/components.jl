# Component-level unit tests of the Stage A Dyad components, one testset per component, against the
# fixtures in dyad/test_harnesses.dyad.
#
# Why these tests are trustworthy (design rules):
#   - Every expected value is a hand calculation from parameters read back off the solution (sol.ps),
#     derived before the simulation ran. `retcode == Success` is never the assertion.
#   - Each fixture isolates one component behind pinned boundary conditions (FixedTemperature, ideal
#     source, lossless capacitor), so a failing check localizes to that component and usually to one
#     term of it. The system-level tests (test/physics.jl) exercise the same components inside the
#     rigs, where a compensating pair of errors — e.g. a doubled R_on against a halved R_load — can
#     cancel; in these fixtures nothing is available to cancel against.
using PowerElectronics, DyadInterface, ModelingToolkit, Test

@testset "ConductionMOSFET: channel law with α₂ and ΔR active" begin
    # Aim: verify every term of R_on = R₀(1 + α₁ΔT_j + α₂ΔT_j²) + ΔR, the off-state branch
    # i = G_off·v, the loss P_loss = v·i, and the heat-port sign (Q_flow = −P_loss: heat LEAVES the
    # die). The fixture pins T_j with FixedTemperature (no electro-thermal feedback) and puts the
    # ideal supply straight across the channel (no load), so v = V_sup in both gate states and the
    # measured current is the channel law alone. α₂ and ΔR are nonzero only here: the rigs leave them
    # at 0, so without this testset a sign/scaling error in either would first surface mid-Stage-A fit.
    # Interpretation: each check isolates one term — currents pass but heat fails → Q_flow sign broke;
    # on-current off by the ratio of R_on values → one of the R terms; off-current wrong → off branch.
    res = MOSFETConduction()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    R₀ = sol.ps[m.mosfet.R₀]; α₁ = sol.ps[m.mosfet.α₁]; α₂ = sol.ps[m.mosfet.α₂]
    ΔR = sol.ps[m.mosfet.ΔR]; G_off = sol.ps[m.mosfet.G_off]
    V = sol.ps[m.V_sup]; t_on = sol.ps[m.t_on]
    ΔT_j = sol.ps[m.T_j_fix] - sol.ps[m.mosfet.T_ref]
    # the test only makes sense if the quadratic and aging terms actually contribute
    @test α₂ * ΔT_j^2 > 0.1 && ΔR > 0
    R_exp = R₀ * (1 + α₁ * ΔT_j + α₂ * ΔT_j^2) + ΔR                    # 0.4766 Ω at the defaults
    t_off_s, t_on_s = t_on / 2, (t_on + 1.0) / 2                       # mid off- and on-intervals
    # off state: exact leakage law, not merely "small" (a wrong conductance would also look small)
    @test sol(t_off_s, idxs = m.mosfet.i) ≈ G_off * V rtol = 1e-9
    @test sol(t_off_s, idxs = m.mosfet.v) ≈ V rtol = 1e-9
    # on state: current and R_on against the hand value, including α₂ and ΔR
    @test sol(t_on_s, idxs = m.mosfet.R_on) ≈ R_exp rtol = 1e-9
    @test sol(t_on_s, idxs = m.mosfet.i) ≈ V / R_exp rtol = 1e-6       # gate settled to 10·(1−e⁻⁵⁰)
    # loss and heat-port sign, in both states: P_loss = v·i and Q_flow = −P_loss
    for t in (t_off_s, t_on_s)
        v, i = sol(t, idxs = m.mosfet.v), sol(t, idxs = m.mosfet.i)
        @test sol(t, idxs = m.mosfet.P_loss) ≈ v * i rtol = 1e-9
        @test sol(t, idxs = m.mosfet.junction.Q_flow) ≈ -v * i rtol = 1e-9
    end
end

@testset "CauerNetwork: heat step, aging scaling, first law" begin
    # Aim: with constant P into port_j and port_c pinned, the ladder's steady state is pure Kirchhoff:
    # T₁−T_sink = P((1+a)R_th[1] + R_th[2]) and T₂−T_sink = P·R_th[2]. Three independent failure modes are
    # covered: (i) initialization — both nodes must start at T₀, the contract Stage A's warm-start
    # uses; (ii) the steady drops — a wrong ladder wiring (resistor between the wrong nodes) moves
    # them; (iii) the integral balance ∫Q_j + ∫Q_c = ΣC_th(T−T₀) — the first law across the component
    # boundary including BOTH port sign conventions. (iii) is the check that catches a lost or doubled
    # flow term, which (ii) cannot see: steady state constrains flow ratios, not their bookkeeping.
    # The a = 0.5 rerun must raise the drop across stage 1 by exactly 1.5 and leave the drop across
    # stage 2 as it is: the aging state acts on the first stage only, R_eff[1] = R_th[1](1+a).
    res = CauerStep()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    P = sol.ps[m.P_in]; T_s = sol.ps[m.T_sink]
    R₁, R₂ = sol.ps[m.jc.R_th[1]], sol.ps[m.jc.R_th[2]]
    C₁, C₂ = sol.ps[m.jc.C_th[1]], sol.ps[m.jc.C_th[2]]
    # (i) initialization: every node starts at T₀
    @test sol(0.0, idxs = m.jc.T[1]) ≈ T_s rtol = 1e-12
    @test sol(0.0, idxs = m.jc.T[2]) ≈ T_s rtol = 1e-12
    # (ii) steady state after ≈50 slow time constants: each stage drop is P·R_th[k]
    T₁e, T₂e = sol[m.jc.T[1]][end], sol[m.jc.T[2]][end]
    @test T₁e ≈ T_s + P * (R₁ + R₂) rtol = 1e-8
    @test T₂e ≈ T_s + P * R₂ rtol = 1e-8
    # (iii) first law: heat in through port_j plus heat in through port_c (negative: it leaves) equals
    # the stored energy. Trapezoid on a dense grid; the fast die mode (τ₁ ≈ 5 ms) sets the grid.
    ts = range(0, sol.t[end], length = 300_001); dt = step(ts)
    Qj = sol(ts, idxs = m.jc.port_j.Q_flow).u; Qc = sol(ts, idxs = m.jc.port_c.Q_flow).u
    trap(v) = dt * (sum(v) - (v[1] + v[end]) / 2)
    @test trap(Qj) + trap(Qc) ≈ C₁ * (T₁e - T_s) + C₂ * (T₂e - T_s) rtol = 1e-3
    # aging on stage 1 only, a = 0.5: R_eff = [1.0·1.5, 2.1] K/W at P = 2 W, so by hand
    # T₁ − T_sink = 2·(1.5 + 2.1) = 7.2 K and T₂ − T_sink = 2·2.1 = 4.2 K, the a = 0 value of node 2.
    # With stage 2 scaled as well the rises would be 2·1.5·3.1 = 9.3 K and 2·1.5·2.1 = 6.3 K.
    res_a = CauerStep(model = PowerElectronics.CauerStepTest(; name = :CauerStepTest, a_net = 0.5))
    sol_a, m_a = res_a.sol, symbolic_container(res_a)
    @test sol_a[m_a.jc.T[1]][end] ≈ T_s + P * (1.5 * R₁ + R₂) rtol = 1e-8
    @test sol_a[m_a.jc.T[2]][end] ≈ T_s + P * R₂ rtol = 1e-8
end

@testset "DutyAveragedHeat: the cycle-averaging contract" begin
    # Aim: the operator's whole contract is port_network.Q_flow = −duty·port_device.Q_flow with
    # port_device.T = port_network.T. Feeding its output into a lossless capacitor turns that into a
    # straight line, T(t) = T_start + duty·P/C·t: no other admissible slope exists, so a full-P leak
    # (missing duty), a sign flip, or a swapped port shows up as the wrong slope or a cooling trace.
    # Interpretation: this component is DELIBERATELY non-conservative — (1−duty)·P must vanish, that
    # is the cycle average it implements — so the right check is the operator contract, not port-to-
    # port energy conservation (which would fail by design). The duty = 0 rerun must hold the
    # capacitor at T_start exactly: that is the path the rig's hysteresis controller uses to
    # switch the heat off, so an off-state leak here becomes a wrong hold temperature there.
    res = DutyAveraging()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    P = sol.ps[m.P_dev]; d = sol.ps[m.duty_k]; C = sol.ps[m.C_sink]; T₀ = sol.ps[m.T_start]
    for t in (1.0, 3.0, 5.0)
        @test sol(t, idxs = m.store.T) ≈ T₀ + d * P / C * t rtol = 1e-8        # the line, no curvature
        @test sol(t, idxs = m.averaging.port_network.Q_flow) ≈ -d * P rtol = 1e-8  # duty × device heat
        @test sol(t, idxs = m.averaging.port_device.Q_flow) ≈ P rtol = 1e-8
        @test sol(t, idxs = m.averaging.port_device.T) ≈ sol(t, idxs = m.store.T) rtol = 1e-12
    end
    # duty = 0: heat fully blocked, capacitor pinned at T_start — the controller-off path
    res₀ = DutyAveraging(model = PowerElectronics.DutyAveragingTest(; name = :DutyAveragingTest, duty_k = 0.0))
    sol₀, m₀ = res₀.sol, symbolic_container(res₀)
    @test maximum(abs.(sol₀[m₀.store.T] .- sol₀.ps[m₀.T_start])) < 1e-9
end

@testset "PWMSource: levels, phase stability, duty fraction" begin
    # Aim: the gate command's shape (offset / offset+amplitude), its timing over many periods, and the
    # duty FRACTION itself, which is the number the thermal model consumes (heat scales with duty).
    # Levels are sampled mid-high and mid-low in periods 1, 5 and 10: a wrong floor() construction
    # accumulates phase drift that grows with period count, which a single-period check cannot see.
    # The integrator turns duty into one number, ∫y dt = offset·t_stop + amplitude·duty·(t_stop −
    # t_start), exact over the 10 whole periods of the window — a few-percent duty error would pass
    # every level check yet bias every Stage A heat input by the same factor. The fixture's parameters
    # are deliberately asymmetric (offset 0.5, start_time 0.05): the rig defaults (0, 0) would mask a
    # source that ignores either. duty = 0 and duty = 1 pin the boundary arithmetic of `frac < duty`,
    # where < vs <= flips the output.
    res = PWMWindow()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    A = sol.ps[m.pwm.amplitude]; f = sol.ps[m.pwm.frequency]; d = sol.ps[m.pwm.duty]
    off = sol.ps[m.pwm.offset]; t₀ = sol.ps[m.pwm.start_time]
    T_pwm = 1 / f; t_stop = sol.t[end]
    @test sol(t₀ / 2, idxs = m.pwm.y) ≈ off rtol = 1e-9                 # before start_time: offset only
    for k in (0, 4, 9)                                                   # periods 1, 5, 10: no phase drift
        @test sol(t₀ + k * T_pwm + d * T_pwm / 2, idxs = m.pwm.y) ≈ off + A rtol = 1e-9   # mid high
        @test sol(t₀ + k * T_pwm + (1 + d) * T_pwm / 2, idxs = m.pwm.y) ≈ off rtol = 1e-9 # mid low
    end
    # duty fraction via the integral over 10 whole periods (solver smears each edge by ≤ dtmax)
    @test sol(t_stop, idxs = m.intg.y) ≈ off * t_stop + A * d * (t_stop - t₀) rtol = 1e-3
    # boundary duties: output must be constant-low / constant-high, integral exact
    for (dv, I_exp) in ((0.0, off * t_stop), (1.0, off * t_stop + A * (t_stop - t₀)))
        res_d = PWMWindow(model = PowerElectronics.PWMSourceTest(; name = :PWMSourceTest, pwm_duty = dv))
        sol_d, m_d = res_d.sol, symbolic_container(res_d)
        @test sol_d(t_stop, idxs = m_d.intg.y) ≈ I_exp rtol = 1e-3
        @test sol_d(0.5, idxs = m_d.pwm.y) ≈ off + dv * A rtol = 1e-9
    end
end
