# Observation-model blocks: CurrentOffset, VdsStepArtefact,
# ThermocoupleLag (dyad/observation.dyad) against their fixtures in dyad/test_harnesses.dyad, and the noise
# helper add_measurement_noise (src/observation.jl).
#
# Design rules (test/components.jl lines 4–11): every expected value is a hand calculation from parameters read
# back off the solution (sol.ps), and `retcode == Success` is never the assertion. Tolerances (logged in the @info
# lines): 1e-6 relative for step responses and window means, 1e-4 relative for sine amplitude and phase, the
# noise helper within 3 standard errors at n = 1e5. Each testset logs its measured worst error with @info.
#
# Fixtures that contain a discontinuous input (step, ramp corner, pulse edge) are built from Julia so the edges
# go into `tstops` (a Dyad TransientAnalysis takes no tstops); the smooth ones run through their analyses.
using PowerElectronics, DyadInterface, ModelingToolkit, OrdinaryDiffEqTsit5, Test
using DyadInterface: symbolic_container, simplify_model
using Statistics: mean, std

"Simplified system of a fixture with top-level parameter keywords (Dyad 3.5: set top-level names only)."
obs_build(fixture; kw...) = simplify_model(fixture(; name = :fix, kw...))

"Solve a fixture from t = 0 to `stop` with Tsit5 at tight tolerances and the given discontinuity times."
obs_solve(sys, stop, tstops) =
    solve(ODEProblem(sys, [], (0.0, stop)), Tsit5(); tstops = tstops, abstol = 1e-13, reltol = 1e-12)

"Composite-Simpson mean of observable `idx` over [t1, t2] from the dense solution (n even)."
function obs_window_mean(sol, idx, t1, t2; n = 2000)
    ts = collect(range(t1, t2, length = n + 1))
    v = sol(ts; idxs = idx).u
    w = [i == 1 || i == n + 1 ? 1 : (iseven(i) ? 4 : 2) for i in 1:n+1]
    sum(w .* v) * (ts[2] - ts[1]) / 3 / (t2 - t1)
end

@testset "Observation model" begin

@testset "CurrentOffset: I_d,meas = I_d + offset_I" begin
    # Aim: the block adds a constant, independent of the current level, and nothing else. A ramp from 0 A to
    # 15 A (about the rig's cold on-state current) through the block: y − u must equal offset_I at every
    # instant to 1e-12 A, at the transparent default 0 and at 0.3 A (inside the 0.2–0.4 A seen in the data). A
    # gain error would make the difference drift along the ramp. The integrator gives the integral identity
    # ∫y dt = ∫u dt + offset_I·t, with ∫u over the ramp a hand calculation.
    worst = 0.0
    for off in (0.0, 0.3)
        sys = obs_build(PowerElectronics.CurrentOffsetTest; off_I = off)
        t0r, dur = 0.1, 0.8
        sol = obs_solve(sys, 1.0, [t0r, t0r + dur])
        @test sol.ps[sys.meas.offset_I] == off                  # the top-level value reached the block
        H = sol.ps[sys.ramp.height]
        @test sol.ps[sys.ramp.start_time] == t0r && sol.ps[sys.ramp.duration] == dur
        ts = range(0.0, 1.0, length = 201)
        d = [sol(t; idxs = sys.meas.y) - sol(t; idxs = sys.ramp.y) for t in ts]
        err = maximum(abs.(d .- off))
        worst = max(worst, err)
        @test err <= 1e-12
        # ∫u dt over [0, 1]: 0 before t0r, the ramp triangle, then the plateau
        ∫u = H * dur / 2 + H * (1.0 - t0r - dur)
        @test sol[sys.intg.y][end] ≈ ∫u + off * 1.0 rtol = 1e-9
    end
    @info "CurrentOffset: max |y − u − offset_I| = $(worst) A (bar 1e-12 A)"
end

@testset "VdsStepArtefact: step response, unit gain, static offset" begin
    # Aim: (i) after a step ΔU at t_step the error is e(t) = k_V·ΔU·exp(−(t − t_step)/τ_V): a wrong jump height
    # points at k_V or the sign of (u − x), a wrong decay rate at τ_V; (ii) at constant input the output equals
    # the input (gain 1) — the fixture starts the filter at the pre-step level, so y = u before the step;
    # the rerun with offset_V = 0.09 V checks y − u − e = offset_V. τ_V = 70 µs ≠ 1, so a missing
    # division by τ_V cannot pass. Tolerance: 1e-6 relative.
    sys = obs_build(PowerElectronics.VdsArtefactStepTest)
    t_s = 2e-4
    sol = obs_solve(sys, 1e-3, [t_s])
    k = sol.ps[sys.art.k_V]; τ = sol.ps[sys.art.τ_V]; ΔU = sol.ps[sys.dU]
    @test sol.ps[sys.t_step] == t_s && sol.ps[sys.art.x_start] == sol.ps[sys.U0]
    @test k > 0 && τ == 70e-6                                   # the block default τ_V
    pre = [abs(sol(t; idxs = sys.art.y) - sol(t; idxs = sys.art.u)) for t in range(0.0, 0.99t_s, length = 100)]
    @test maximum(pre) <= 1e-12                                 # (ii) gain 1 at constant input
    rels = Float64[]
    for s in (0.01, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0) .* τ
        e = sol(t_s + s; idxs = sys.art.e)
        push!(rels, abs(e / (k * ΔU * exp(-s / τ)) - 1))
    end
    @test maximum(rels) <= 1e-6                                 # (i) jump and decay
    # the output is input + offset + error, with offset 0 here
    @test sol(t_s + τ; idxs = sys.art.y) ≈ sol(t_s + τ; idxs = sys.art.u) + sol(t_s + τ; idxs = sys.art.e) rtol = 1e-12
    # static offset: 0.09 V between the voltage channels
    sys9 = obs_build(PowerElectronics.VdsArtefactStepTest; off_V = 0.09)
    sol9 = obs_solve(sys9, 1e-3, [t_s])
    @test sol9.ps[sys9.art.offset_V] == 0.09
    off_err = maximum(abs(sol9(t; idxs = sys9.art.y) - sol9(t; idxs = sys9.art.u) - sol9(t; idxs = sys9.art.e) - 0.09)
                      for t in range(0.0, 1e-3, length = 101))
    @test off_err <= 1e-12
    @test maximum(abs(sol9(t; idxs = sys9.art.y) - sol9(t; idxs = sys9.art.u) - 0.09)
                  for t in range(0.0, 0.99t_s, length = 50)) <= 1e-12
    @info "VdsStepArtefact step: max rel error of e(t) = $(maximum(rels)) (bar 1e-6); max |y − u| before the step = $(maximum(pre)) V; offset check $(off_err) V"
end

@testset "VdsStepArtefact: windows of the 1 kHz, 40 % gate drive" begin
    # Aim: on the NASA gate drive the mean error over the export's two on-state windows,
    # vds_early (20–30 % of the on-run) and vds_on (78–98 %), equals the
    # closed-form periodic steady state: z₀ = a(1 − e^{−T_off/τ})/(1 − e^{−T/τ}) just after turn-on (a the
    # turn-on step), and mean e over [f₁, f₂]·T_on = k·z₀·τ·(e^{−f₁T_on/τ} − e^{−f₂T_on/τ})/((f₂ − f₁)T_on).
    # This is the quantity compared with the data. Period 6 is used: the start transient decays
    # by e^{−T/τ} ≈ 6e-7 per period. Tolerance: 1e-6 relative.
    sys = obs_build(PowerElectronics.VdsArtefactPWMTest)
    f = 1000.0; dty = 0.4                                       # fixture values, checked against sol.ps below
    T = 1 / f
    edges = sort(vcat([j * T for j in 1:5], [j * T + dty * T for j in 0:5]))
    sol = obs_solve(sys, 6T, edges)
    @test sol.ps[sys.gate.frequency] == f && sol.ps[sys.gate.duty] == dty
    k = sol.ps[sys.art.k_V]; τ = sol.ps[sys.art.τ_V]
    a = sol.ps[sys.V_on] - sol.ps[sys.V_off]
    @test sol.ps[sys.gate.amplitude] ≈ a rtol = 1e-15
    Ton = dty * T; Toff = T - Ton
    z0 = a * (1 - exp(-Toff / τ)) / (1 - exp(-T / τ))
    closed(f1, f2) = k * z0 * τ * (exp(-f1 * Ton / τ) - exp(-f2 * Ton / τ)) / ((f2 - f1) * Ton)
    t0 = 5T
    # the pulse levels in period 6 (the input is what the closed form assumes)
    @test sol(t0 + 0.5Ton; idxs = sys.art.u) == sol.ps[sys.V_on]
    @test sol(t0 + Ton + 0.5Toff; idxs = sys.art.u) == sol.ps[sys.V_off]
    rels = Dict{String, Float64}()
    for (nm, f1, f2) in (("vds_early", 0.20, 0.30), ("vds_on", 0.78, 0.98))
        m = obs_window_mean(sol, sys.art.e, t0 + f1 * Ton, t0 + f2 * Ton)
        rels[nm] = abs(m / closed(f1, f2) - 1)
        @test rels[nm] <= 1e-6
    end
    @info "VdsStepArtefact PWM windows: rel error vds_early = $(rels["vds_early"]), vds_on = $(rels["vds_on"]) (bar 1e-6); closed-form means $(closed(0.2, 0.3)) V, $(closed(0.78, 0.98)) V at k_V = $k"
end

@testset "ThermocoupleLag: step response" begin
    # Aim: T_meas = T0 before the step and T0 + dT·(1 − e^{−(t − t_step)/τ_tc}) after it, from the start
    # parameter T_start = T0. τ_tc = 0.37 s ≠ 1, so a missing division by τ_tc cannot pass.
    # Tolerance: 1e-6 relative.
    sys = obs_build(PowerElectronics.LagStepTest)
    t_s = 1.0
    sol = obs_solve(sys, 5.0, [t_s])
    τ = sol.ps[sys.lag.τ_tc]; T0 = sol.ps[sys.T0]; dT = sol.ps[sys.dT]
    @test sol.ps[sys.t_step] == t_s && τ == sol.ps[sys.tau] && τ != 1.0
    @test sol(0.0; idxs = sys.lag.T_meas) == sol.ps[sys.lag.T_start] == T0     # starts at the parameter
    pre = maximum(abs(sol(t; idxs = sys.lag.y) - T0) for t in range(0.0, 0.99t_s, length = 100))
    @test pre <= 1e-9
    rels = [abs((sol(t_s + s; idxs = sys.lag.y) - T0) / (dT * (1 - exp(-s / τ))) - 1) for s in (0.05, 0.25, 0.5, 1, 2, 4, 8) .* τ]
    @test maximum(rels) <= 1e-6
    @info "ThermocoupleLag step: max rel error of the step response = $(maximum(rels)) (bar 1e-6)"
end

@testset "ThermocoupleLag: sine at 4.85 s, amplitude, phase, mean" begin
    # Aim: the first-order lag's frequency response at test 30's measured burst period: amplitude ratio
    # 1/√(1 + (ωτ)²) and phase lag atan(ωτ) (tolerance 1e-4 relative), and unit gain at zero frequency: the mean
    # of T_meas over whole periods equals the input mean (1e-6 relative), the property a cycle-averaged model
    # relies on. Read over periods 11–20 (transient below e^{−48}) by projection on sin and cos; the periodic
    # trapezoid rule is exact to rounding for a band-limited periodic integrand.
    r = PowerElectronics.LagSine()
    sol, m = r.sol, symbolic_container(r)
    τ = sol.ps[m.lag.τ_tc]; A = sol.ps[m.A]; Tm = sol.ps[m.T_mean]; P = sol.ps[m.period]
    @test sol.ps[m.src.frequency] ≈ 1 / P rtol = 1e-15
    ω = 2π / P; N = 4000
    ts = collect(range(10P, 20P, length = 10N + 1)[1:end-1])
    v = sol(ts; idxs = m.lag.T_meas).u
    u = sol(ts; idxs = m.src.y).u
    μ = mean(v)
    a = 2 * mean((v .- μ) .* sin.(ω .* ts)); b = 2 * mean((v .- μ) .* cos.(ω .* ts))
    ratio = hypot(a, b) / A; φ = atan(-b, a)
    rel_amp = abs(ratio * sqrt(1 + (ω * τ)^2) - 1)
    rel_phase = abs(φ / atan(ω * τ) - 1)
    @test rel_amp <= 1e-4
    @test rel_phase <= 1e-4
    @test abs(mean(u) - Tm) <= 1e-9 * Tm                        # the input itself has mean T_mean over the window
    rel_mean = abs(μ - Tm) / Tm
    @test rel_mean <= 1e-6
    @test abs(μ - Tm) <= 1e-6 * A                               # also tight against the amplitude
    @info "ThermocoupleLag sine: rel error amplitude = $(rel_amp), phase = $(rel_phase) (bar 1e-4); mean rel error = $(rel_mean) (bar 1e-6), |Δmean| = $(abs(μ - Tm)) K"
end

@testset "ThermocoupleLag: attached block leaves the network unchanged" begin
    # Aim: the block reads T through a TemperatureSensor (port.Q_flow = 0), so the thermal states are identical
    # with and without it. Twin nodes in one solve — a read through sensor + lag, b bare — must coincide, and
    # both must follow T_amb + ΔT·e^{−t/τ_th}. The lag's reading of that exponential is the hand calculation
    # T_amb + ΔT·(τ_th e^{−t/τ_th} − τ_tc e^{−t/τ_tc})/(τ_th − τ_tc). Tolerance: 1e-6 relative.
    r = PowerElectronics.LagTwin()
    sol, m = r.sol, symbolic_container(r)
    C = sol.ps[m.C_node]; R = sol.ps[m.R_node]; Th = sol.ps[m.T_hot]; Ta = sol.ps[m.T_amb]; τc = sol.ps[m.lag.τ_tc]
    τh = R * C; ΔT = Th - Ta
    @test sol.ps[m.node_a.C] == sol.ps[m.node_b.C] == C && sol.ps[m.r_a.R] == sol.ps[m.r_b.R] == R
    same = maximum(abs.(sol[m.node_a.T] .- sol[m.node_b.T]) ./ Th)
    @test same <= 1e-12
    @test maximum(abs.(sol[m.sensor.port.Q_flow])) == 0.0
    ts = range(0.0, sol.t[end], length = 601)
    rel_a = maximum(abs((sol(t; idxs = m.node_a.T) - Ta) / (ΔT * exp(-t / τh)) - 1) for t in ts)
    rel_b = maximum(abs((sol(t; idxs = m.node_b.T) - Ta) / (ΔT * exp(-t / τh)) - 1) for t in ts)
    Tmeas(t) = Ta + ΔT * (τh * exp(-t / τh) - τc * exp(-t / τc)) / (τh - τc)
    rel_m = maximum(abs((sol(t; idxs = m.lag.T_meas) - Ta) / (Tmeas(t) - Ta) - 1) for t in ts)
    @test rel_a <= 1e-6 && rel_b <= 1e-6
    @test rel_m <= 1e-6
    @info "ThermocoupleLag twin: max rel |T_a − T_b| = $(same); rel error of the rises a, b = $(rel_a), $(rel_b); of the reading = $(rel_m) (bar 1e-6)"
end

@testset "add_measurement_noise: seeded Gaussian noise" begin
    # Aim: the helper adds N(0, σ²) noise to saved samples, reproducibly. At n = 1e5 the sample mean and sd of
    # the added noise lie within 3 standard errors of 0 and σ (SE σ/√n and ≈ σ/√(2(n − 1))). The same seed gives
    # the same series, another seed another; the input is not modified; σ = 0 returns the input; a negative σ
    # is refused. σ = 0.0036 K is the pooled σ_Tc, the case-temperature noise estimated from the data.
    n = 10^5; σ = 0.0036; seed = 20261031
    x = collect(range(483.0, 484.0, length = n))
    x0 = copy(x)
    y = PowerElectronics.add_measurement_noise(x, σ; seed = seed)
    ε = y .- x
    z_mean = abs(mean(ε)) / (σ / sqrt(n))
    z_sd = abs(std(ε) - σ) / (σ / sqrt(2 * (n - 1)))
    @test z_mean <= 3
    @test z_sd <= 3
    @test y == PowerElectronics.add_measurement_noise(x, σ; seed = seed)
    @test y != PowerElectronics.add_measurement_noise(x, σ; seed = seed + 1)
    @test x == x0
    @test PowerElectronics.add_measurement_noise(x, 0.0; seed = seed) == x
    @test_throws ArgumentError PowerElectronics.add_measurement_noise(x, -1.0; seed = seed)
    @info "noise helper: |mean|/SE = $(z_mean), |sd − σ|/SE = $(z_sd) at n = $n (bar 3)"
end

end
