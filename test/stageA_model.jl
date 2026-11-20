# Stage A model: junction steady state and aging-law placement, and the energy balance over
# test 30's s1. Measured headroom at the time of writing: energy imbalance ~2e-11 against 1e-4;
# junction steady state ~3e-16 against 1e-6.
using Test, CSV, DataFrames
using ModelingToolkit, OrdinaryDiffEqDefault
using DyadInterface: symbolic_container, simplify_model
using DyadData: DyadTimeseries

@testset "Stage A model" begin
    # test 30's 2N-fixture case-side parameters, as in StageATest30
    C30, Gc30, Rcf30, Cf30 = 1.9817297999882928, 0.009401386147471817, 42.21091872296592, 1.1088981194041447
    T0K = 298.596                # s1's first measured T_c, = T_start = T_amb of StageATest30
    stop = 1741.008              # s1 end, = stop of StageATest30

    @testset "junction steady state and aging law" begin
        r = JunctionStep()
        sol = r.sol
        m = symbolic_container(r)
        @test string(sol.retcode) == "Success"
        # steady state is the hand calculation T_case + P_in·R_jc·(1 + a)
        @test abs(sol[m.jc.T_j][end] - 406.2) <= 1e-6 * 406.2
        # a = 0.5 rebuilt from Julia: the rise must scale by exactly 1.5 — the aging law
        # sits on R, not on C (the same contract CauerStep pins for CauerNetwork)
        @named fix = PowerElectronics.JunctionStepTest(a_net = 0.5)
        sys = simplify_model(fix)
        sol2 = solve(ODEProblem(sys, [], (0.0, 0.8)); abstol = 1e-10, reltol = 1e-10)
        @test string(sol2.retcode) == "Success"
        @test abs(sol2[sys.jc.T_j][end] - 409.3) <= 1e-6 * 409.3
    end

    @testset "two-stage junction: the aging state scales stage 1 only" begin
        # StageAJunction2 with constant heat P = 2 W and the case held at 400 K (two-row series),
        # R_jc1 = 1.0 K/W, R_jc2 = 2.1 K/W, a = 0.5. By hand, with a on stage 1 only:
        # die = 400 + P·(R_jc1·(1 + a) + R_jc2) = 400 + 2·(1.0·1.5 + 2.1) = 407.2 K and
        # stage 2 = 400 + P·R_jc2 = 400 + 2·2.1 = 404.2 K. With stage 2 scaled as well the two
        # would be 409.3 K and 406.3 K; with no stage scaled, 406.2 K and 404.2 K.
        # 10 s ≈ 40 slow time constants (≈ R_jc2·(C_die1 + C_die2) = 0.23 s).
        heat = DyadTimeseries([0.0 2.0; 10.0 2.0]; independent_var = "time", dependent_vars = ["P"])
        tc = DyadTimeseries([0.0 400.0; 10.0 400.0]; independent_var = "time", dependent_vars = ["T_c"])
        @named two = PowerElectronics.StageAJunction2(heat_data = heat, tc_data = tc, T_start = 400.0,
                                                      R_jc1 = 1.0, R_jc2 = 2.1, a = 0.5)
        sys2 = simplify_model(two)
        sol3 = solve(ODEProblem(sys2, [], (0.0, 10.0)); abstol = 1e-10, reltol = 1e-10)
        @test string(sol3.retcode) == "Success"
        @test abs(sol3[sys2.T_j][end] - 407.2) <= 1e-6 * 407.2
        @test abs(sol3[sys2.jc2.T_j][end] - 404.2) <= 1e-6 * 404.2
    end

    @testset "energy balance, test 30 s1" begin
        r = StageATest30()
        sol = r.sol
        m = symbolic_container(r)
        @test string(sol.retcode) == "Success"

        # energy balance: ∫P dt (analytic from the step series, clipped to stop) equals the
        # stored-heat change of the three nodes plus the integrated ambient-path heat E_amb
        h = CSV.read(joinpath(pkgdir(PowerElectronics), "assets", "stageA", "test30_heat.csv"), DataFrame)
        E_in = sum(h.P[i] * (min(h.time[i+1], stop) - h.time[i]) for i in 1:nrow(h)-1 if h.time[i] < stop)
        ΔE = C30 * (sol[m.T_c][end] - T0K) + 1e-4 * (sol[m.T_j][end] - T0K) + Cf30 * (sol[m.T_f][end] - T0K)
        @test abs(E_in - ΔE - sol[m.E_amb][end]) <= 1e-4 * E_in
    end
end
