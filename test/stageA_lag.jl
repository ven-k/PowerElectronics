# StageAModelLag (dyad/observation.dyad): StageAModel with ThermocoupleLag on its T_c output.
#
# Design rules (test/components.jl lines 4–11): every expected value is a hand calculation from parameters read back
# off the solution (sol.ps), and `retcode == Success` is never the assertion. Each testset logs its worst error with
# @info. Systems are built from Julia, with top-level parameter keywords only (Dyad 3.5),
# driven by test 36's measured burst heat, and every solve passes the burst edges as tstops.
#
# Checks: (1) the lag is a signal block and draws no heat, so T_j, T_c, T_f equal a plain StageAModel solve with
# the same parameters, same solver, tolerances and tstops (≤ 1e-9 relative); (2) the reading starts at the
# parameter T_start; (3) the reading obeys dT_c,meas/dt = (T_c − T_c,meas)/tau_tc: over a window,
# T_c,meas(t2) − T_c,meas(t1) = ∫ (T_c − T_c,meas)/tau_tc dt, the integral by composite Simpson on the dense
# solution split at the burst edges (≤ 1e-6 relative). tau_tc = 2.5 s, so a missing division (τ = 1) fails by 150 %.
using PowerElectronics, DyadInterface, ModelingToolkit, OrdinaryDiffEqTsit5, DyadData, Test
using DyadInterface: simplify_model

"Composite Simpson of f over [a, b], split at the breakpoints inside (n even, per segment)."
function lag_simpson(f, a, b, breaks; n = 2000)
    edges = vcat(a, filter(x -> a < x < b, breaks), b)
    s = 0.0
    for k in 1:length(edges)-1
        x = range(edges[k], edges[k+1], length = n + 1)
        v = f.(x)
        s += step(x) / 3 * (v[1] + v[end] + 4 * sum(v[2:2:end-1]) + 2 * sum(v[3:2:end-2]))
    end
    s
end

@testset "StageAModelLag" begin
    heat = DyadData.DyadTimeseries("dyad://PowerElectronics/stageA/test36_heat.csv";
                                   independent_var = "time", dependent_vars = ["P"])
    # Non-default values, near test 36's fitted values (see Rig36), so a dropped keyword shows in the read-back.
    kw = (T_start = 301.0, T_amb = 300.3, C = 0.9, Gc = 0.0080, R_cf = 12.5, C_f = 1.6, R_jc = 4.7, C_die = 0.063)
    tau = 2.5
    sysL = simplify_model(PowerElectronics.StageAModelLag(; name = :lagm, heat_data = heat, tau_tc = tau, kw...))
    sys0 = simplify_model(PowerElectronics.StageAModel(; name = :plain, heat_data = heat, kw...))
    T_END = 200.0
    ts = filter(t -> 0 < t < T_END, PowerElectronics.stageA_tstops(36))
    lag_solve(sys) = solve(ODEProblem(sys, [], (0.0, T_END)), Tsit5(); tstops = ts, abstol = 1e-12, reltol = 1e-12)
    solL, sol0 = lag_solve(sysL), lag_solve(sys0)

    @testset "parameters read back" begin
        # Aim: both solves carry the values set through the top-level names, and the lag block's own parameters
        # are bound to them (lag.τ_tc = tau_tc, lag.T_start = T_start).
        for (k, v) in pairs(kw)
            @test solL.ps[getproperty(sysL, k)] == v
            @test sol0.ps[getproperty(sys0, k)] == v
        end
        @test solL.ps[sysL.tau_tc] == tau
        @test solL.ps[sysL.lag.τ_tc] == solL.ps[sysL.tau_tc]
        @test solL.ps[sysL.lag.T_start] == solL.ps[sysL.T_start]
        @test abs(solL.ps[sysL.tau_tc] - 1) >= 1      # a missing division by tau_tc is visible below
        @info "StageAModelLag read-back: tau_tc = $(solL.ps[sysL.tau_tc]) s, lag.τ_tc = $(solL.ps[sysL.lag.τ_tc]) s"
    end

    @testset "lag draws no heat: T_j, T_c, T_f equal StageAModel" begin
        # Aim: attaching the signal block leaves the thermal network unchanged; the reference is the plain model at
        # the same parameters. Bar 1e-9 relative (both solves at 1e-12).
        tt = range(0.0, T_END, length = 4001)
        worst = 0.0
        for (a, b) in ((sysL.T_j, sys0.T_j), (sysL.T_c, sys0.T_c), (sysL.T_f, sys0.T_f))
            ya, yb = solL(tt; idxs = a).u, sol0(tt; idxs = b).u
            worst = max(worst, maximum(abs.(ya .- yb) ./ abs.(yb)))
        end
        @test worst <= 1e-9
        # the comparison is not trivial: the network heats by > 100 K and the reading lags T_c by > 1 K
        rise = maximum(sol0[sys0.T_c]) - sol0.ps[sys0.T_start]
        gap = maximum(abs.(solL(tt; idxs = sysL.T_c).u .- solL(tt; idxs = sysL.T_c_meas).u))
        @test rise > 100 && gap > 1
        @info "StageAModelLag vs StageAModel: max rel |ΔT| over T_j, T_c, T_f = $(worst) (bar 1e-9); T_c rise $(rise) K, max |T_c − T_c_meas| $(gap) K"
    end

    @testset "reading starts at T_start" begin
        rel = abs(solL[sysL.T_c_meas][1] - solL.ps[sysL.T_start]) / solL.ps[sysL.T_start]
        @test rel <= 1e-12
        @info "StageAModelLag T_c_meas(0) = $(solL[sysL.T_c_meas][1]) K vs T_start = $(solL.ps[sysL.T_start]) K (rel $(rel))"
    end

    @testset "reading obeys dT_c,meas/dt = (T_c − T_c,meas)/tau_tc" begin
        # Aim: integral form of the lag ODE over windows that contain burst edges; quadrature split at the edges
        # (T_c − T_c,meas is C¹ with jumps in its second derivative there). Bar 1e-6 relative.
        τ = solL.ps[sysL.tau_tc]
        worst = 0.0
        for (t1, t2) in ((20.0, 60.0), (100.0, 140.0))
            lhs = solL(t2; idxs = sysL.T_c_meas) - solL(t1; idxs = sysL.T_c_meas)
            rhs = lag_simpson(x -> (solL(x; idxs = sysL.T_c) - solL(x; idxs = sysL.T_c_meas)) / τ, t1, t2, ts)
            worst = max(worst, abs(lhs - rhs) / abs(lhs))
        end
        @test worst <= 1e-6
        @info "StageAModelLag lag ODE, integral form: worst rel error $(worst) over 20–60 s and 100–140 s (bar 1e-6)"
    end
end
