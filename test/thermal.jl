# CaseToAmbient (dyad/thermal.dyad) through its CaseCooldown analysis: a lumped case cooling from 250 °C with no
# heat input. Reads the parameters from the solution, so the checks hold for any positive defaults.
using PowerElectronics, DyadInterface, ModelingToolkit, Test

@testset "CaseToAmbient free cool-down" begin
    res = CaseCooldown()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    C, Ta = sol.ps[m.C_c], sol.ps[m.T_amb]
    Gc, Gr, δ = sol.ps[m.to_amb.Gc], sol.ps[m.to_amb.Gr], sol.ps[m.to_amb.δ]
    ts = range(0, sol.t[end], length = 120_001)
    T, Q = sol(ts, idxs = m.T_c).u, sol(ts, idxs = m.to_amb.Q_flow).u
    # the element implements the stated law: convection Gc ΔT (ΔT² + δ²)^(1/8) plus radiation Gr σ (T⁴ − T∞⁴)
    law(T) = Gc * (T - Ta) * ((T - Ta)^2 + δ^2)^0.125 + Gr * 5.670374419e-8 * (T^4 - Ta^4)
    @test maximum(abs.(Q .- law.(T))) < 1e-8 * maximum(abs.(Q))
    # energy: heat released by the case = heat passed to ambient (trapezoid rule)
    E_out = step(ts) * (sum(Q) - (Q[1] + Q[end]) / 2)
    @test C * (T[1] - T[end]) ≈ E_out rtol = 1e-4
    # nonlinear path: the apparent time constant τ = C (T − T∞)/Q grows as the case cools. A constant R_ca
    # would keep it fixed. Measured on test 36: 58 s at 135 °C, 183 s at 47 °C.
    τ(T) = C * (T - Ta) / law(T)
    @test T[end] < 273.15 + 60                                  # the run reaches the low-temperature end
    @test τ(273.15 + 50) > 1.2 * τ(273.15 + 150) > 0
end

@testset "CaseToAmbient: reversed gradient (case colder than ambient)" begin
    # Aim: the sign behaviour of the δ-smoothed convection law. The element is written as an odd
    # function of ΔT = T_a − T_b; the cool-down test above only ever sees ΔT > 0, so a law that was
    # accidentally even (e.g. |ΔT|^1.25 with a dropped sign) would pass it and then push heat the
    # wrong way the first time the rig model starts below ambient. Here the case starts 40 K BELOW
    # ambient: Q_flow (positive a → b) must be negative throughout, the case must warm monotonically
    # toward ambient and must NOT overshoot — this is a first-order node, so any overshoot means the
    # heat flow crossed zero with the wrong sign near ΔT = 0, exactly where the δ-smoothing acts.
    # The energy check (C·ΔT_case = −∫Q dt) re-verifies the first law on this branch of the law, and
    # the final gap < 0.5 K confirms the smoothed law still drives ΔT → 0 (δ must flatten the
    # derivative at 0, not create a dead band that strands the state short of ambient).
    mod = PowerElectronics.CaseCooldownTest(; name = :CaseCooldownTest, T_start = 283.15, T_amb = 323.15)
    res = CaseCooldown(model = mod)
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    C, Ta = sol.ps[m.C_c], sol.ps[m.T_amb]
    Gc, Gr, δ = sol.ps[m.to_amb.Gc], sol.ps[m.to_amb.Gr], sol.ps[m.to_amb.δ]
    ts = range(0, sol.t[end], length = 120_001)
    T, Q = sol(ts, idxs = m.T_c).u, sol(ts, idxs = m.to_amb.Q_flow).u
    law(T) = Gc * (T - Ta) * ((T - Ta)^2 + δ^2)^0.125 + Gr * 5.670374419e-8 * (T^4 - Ta^4)
    @test maximum(abs.(Q .- law.(T))) < 1e-10                   # the law holds with ΔT < 0 too
    @test all(<(0), Q)                                          # heat flows b → a: INTO the cold case
    @test all(>=(0), diff(T))                                   # monotone warming…
    @test maximum(T) <= Ta                                      # …without overshooting ambient
    E_in = -step(ts) * (sum(Q) - (Q[1] + Q[end]) / 2)           # trapezoid; Q < 0, so −∫Q is heat gained
    @test C * (T[end] - T[1]) ≈ E_in rtol = 1e-4
    @test Ta - T[end] < 0.5                                     # the smoothed law reaches ambient
end

@testset "CaseToAmbient: radiation branch with Gr > 0" begin
    # Aim: the Stefan–Boltzmann term with a NONZERO coefficient. The fitted default is Gr = 0 (not
    # identifiable from the rig data), so every other test in this file exercises Q_rad
    # only at zero — a sign error or a dropped σ in the radiation term would be invisible until
    # someone sets Gr from a datasheet. Gr here is an earlier fit value, i.e. the TO-220
    # geometric scale (emissivity × area × view factor), so the term contributes a physical ~30 % of
    # the initial heat flow rather than a numerically negligible sliver.
    # Interpretation: Q_rad must equal Gr·σ·(T⁴ − T∞⁴) pointwise (the observed variable is the term
    # itself, so this pins sign, σ, and the fourth powers); it must be positive while the case is hot;
    # and the case must cool strictly faster than the Gr = 0 default at the same instant — radiation
    # only ever adds a parallel loss path, it cannot slow the approach to ambient.
    mod = PowerElectronics.CaseCooldownTest(; name = :CaseCooldownTest, Gr_ca = 6.108e-4)
    res = CaseCooldown(model = mod)
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    Gc, Gr, δ = sol.ps[m.to_amb.Gc], sol.ps[m.to_amb.Gr], sol.ps[m.to_amb.δ]
    Ta = sol.ps[m.T_amb]
    @test Gr > 0                                                # the branch under test is actually active
    ts = range(0, sol.t[end], length = 120_001)
    T = sol(ts, idxs = m.T_c).u
    Q_rad = sol(ts, idxs = m.to_amb.Q_rad).u
    Q = sol(ts, idxs = m.to_amb.Q_flow).u
    σSB = 5.670374419e-8
    @test maximum(abs.(Q_rad .- Gr * σSB .* (T .^ 4 .- Ta^4))) < 1e-10
    @test all(>(0), Q_rad[T .> Ta .+ 1.0])                      # radiates outward while hot
    law(T) = Gc * (T - Ta) * ((T - Ta)^2 + δ^2)^0.125 + Gr * σSB * (T^4 - Ta^4)
    @test maximum(abs.(Q .- law.(T))) < 1e-10                   # total still convection + radiation
    # radiation is a parallel loss path: strictly faster cooling than the Gr = 0 default
    res₀ = CaseCooldown()
    sol₀, m₀ = res₀.sol, symbolic_container(res₀)
    @test sol(600.0, idxs = m.T_c) < sol₀(600.0, idxs = m₀.T_c) - 0.5
end
