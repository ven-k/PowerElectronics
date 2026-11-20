# ObservationHarness (dyad/observation.dyad): the junction stage driven by
# test 30's measured s1 T_c, heated by duty_on·I_d²·R_ds during the measured bursts, observed through CurrentOffset
# and VdsStepArtefact (k_V = 0).
#
# Design rules (test/components.jl lines 4–11): every expected value is a hand calculation from parameters read
# back off the solution (sol.ps) and from the data files themselves (the heat step series and the s1 T_c,
# interpolated here by hand), and `retcode == Success` is never the assertion. Each testset logs its worst error
# with @info. Bars, fixed before the first run: 1e-12 relative for algebraic identities at sampled instants
# (they hold to round-off), 1e-9 K for the case temperature (the data's linear interpolant), 1e-6 of the
# window's heat input for the junction energy balance (dense-output quadrature; the window-mean tolerance of test/observation.jl).
#
# The system is built from Julia and solved over the first 400 s of test 30 (bursts and gaps) with
# the burst edges as tstops. Two parameter sets: the defaults, and nonzero offsets together with
# nonzero a, dR and alpha2 so that every parameter's placement is exercised.
using PowerElectronics, DyadInterface, ModelingToolkit, OrdinaryDiffEqDefault, Test
using DyadInterface: simplify_model

@testset "ObservationHarness" begin
    stop = 400.0
    sys = simplify_model(PowerElectronics.ObservationHarness(; name = :oh))
    edges = filter(<=(stop), PowerElectronics.stageA_tstops(30))
    oh_run(op) = solve(ODEProblem(sys, op, (0.0, stop)); tstops = edges, abstol = 1e-10, reltol = 1e-10)
    # test values, not data: offset_I inside the 0.2–0.4 A seen in the data, offset_V the ≈ 0.09 V seen on tests 31–33;
    # a, dR, alpha2 arbitrary nonzero values that make the R_jc·(1 + a), + dR and quadratic terms visible
    sols = (default = oh_run([]),
            shifted = oh_run([sys.offset_I => 0.3, sys.offset_V => 0.09, sys.a => 0.2, sys.dR => 0.01, sys.alpha2 => 2e-6]))

    # the data by hand: heat step series (burst switch) and the s1 case temperature (linear interpolation)
    H = PowerElectronics.DyadTimeseries("dyad://PowerElectronics/stageA/test30_heat.csv"; independent_var = "time", dependent_vars = ["P"])
    S1 = PowerElectronics.DyadTimeseries("dyad://PowerElectronics/stageA/test30_s1.csv"; independent_var = "time", dependent_vars = ["T_c"])
    th, P = Float64.(H["time"]), Float64.(H["P"])
    tc_t, tc = Float64.(S1["time"]), Float64.(S1["T_c"])
    burst(t) = P[searchsortedlast(th, t)] > 1e-9           # step series: value of the last edge at or before t
    function Tc_hand(t)
        i = clamp(searchsortedlast(tc_t, t), 1, length(tc_t) - 1)
        tc[i] + (tc[i+1] - tc[i]) * (t - tc_t[i]) / (tc_t[i+1] - tc_t[i])
    end
    # hand channel law, R0·(1 + α1·ΔT + α2·ΔT²) + dR, from the parameters read off the solution
    pars(sol) = (; (k => sol.ps[getproperty(sys, k)] for k in
                    (:R0, :alpha1, :alpha2, :T_ref, :dR, :R_load, :V_sup, :duty_on, :R_jc, :C_die, :a, :T_start,
                     :offset_I, :offset_V))...)
    Rh(p, Tj) = p.R0 * (1 + p.alpha1 * (Tj - p.T_ref) + p.alpha2 * (Tj - p.T_ref)^2) + p.dR
    # sample instants: a uniform grid that avoids the burst edges, plus the midpoint of every edge interval
    grid = filter(t -> minimum(abs.(edges .- t)) > 1e-6, collect(range(0.05, stop - 0.05, length = 2001)))
    mids = [(edges[i] + edges[i+1]) / 2 for i in 1:length(edges)-1]
    ts = sort(vcat(grid, mids))

    @testset "start and parameters" begin
        worst = 0.0
        for (nm, sol) in pairs(sols)
            p = pars(sol)
            @test sol[sys.T_j][1] == p.T_start                         # T_j(0) = T_start
            @test sol.ps[sys.vds_obs.k_V] == 0.0                       # late-window features: no artefact
            @test sol.ps[sys.id_obs.offset_I] == p.offset_I && sol.ps[sys.vds_obs.offset_V] == p.offset_V
            @test sol.ps[sys.jc.R_jc] == p.R_jc && sol.ps[sys.jc.C_die] == p.C_die && sol.ps[sys.jc.a] == p.a
            # the artefact filter starts at V_ds(0) = V_sup·R(T_start)/(R_load + R(T_start))
            R0h = Rh(p, p.T_start)
            x0 = p.V_sup * R0h / (p.R_load + R0h)
            err = abs(sol[sys.vds_obs.x][1] / x0 - 1)
            worst = max(worst, err)
            @test err <= 1e-12
        end
        @test sols.shifted.ps[sys.offset_I] == 0.3 && sols.shifted.ps[sys.a] == 0.2   # the overrides arrived
        @info "ObservationHarness start: max rel |x(0)/V_ds(0) − 1| = $(worst) (bar 1e-12)"
    end

    @testset "circuit: KVL, channel law, on-state V_ds" begin
        worst = 0.0
        for (nm, sol) in pairs(sols)
            p = pars(sol)
            Tj = sol(ts; idxs = sys.T_j).u
            Id = sol(ts; idxs = sys.I_d).u
            Rd = sol(ts; idxs = sys.R_ds).u
            Vd = sol(ts; idxs = sys.V_ds).u
            Rhand = Rh.(Ref(p), Tj)
            e1 = maximum(abs.(Id .* (p.R_load .+ Rhand) ./ p.V_sup .- 1))   # I_d·(R_load + R_ds) = V_sup
            e2 = maximum(abs.(Rd ./ Rhand .- 1))                              # R_ds = channel law + dR at T_j
            e3 = maximum(abs.(Vd ./ (Id .* Rhand) .- 1))                      # V_ds = I_d·R_ds
            worst = max(worst, e1, e2, e3)
            @test e1 <= 1e-12
            @test e2 <= 1e-12
            @test e3 <= 1e-12
        end
        @info "ObservationHarness circuit: max rel error (KVL, channel law, V_ds) = $(worst) (bar 1e-12)"
    end

    @testset "observation blocks: I_d,meas − I_d = offset_I, V_ds,meas − I_d·R_ds = offset_V" begin
        worst = 0.0
        for (nm, sol) in pairs(sols)
            p = pars(sol)
            Tj = sol(ts; idxs = sys.T_j).u
            Id = sol(ts; idxs = sys.I_d).u
            dI = sol(ts; idxs = sys.I_d_meas).u .- Id
            dV = sol(ts; idxs = sys.V_ds_meas).u .- Id .* Rh.(Ref(p), Tj)
            eI = maximum(abs.(dI .- p.offset_I))
            eV = maximum(abs.(dV .- p.offset_V))
            worst = max(worst, eI, eV)
            @test eI <= 1e-12 * p.V_sup / p.R_load                 # 1e-12 of the largest current, A
            @test eV <= 1e-12 * p.V_sup                             # 1e-12 of the supply, V
        end
        @info "ObservationHarness observation: max |meas − true − offset| = $(worst) (A or V; bar ≈ 1e-11)"
    end

    @testset "bursts: duty_on·I_d²·R_ds while the heat series is positive, zero between" begin
        worst = 0.0
        for (nm, sol) in pairs(sols)
            p = pars(sol)
            on = burst.(mids)
            @test count(on) > 10 && count(.!on) > 10                 # the window holds bursts and gaps
            Tj = sol(mids; idxs = sys.T_j).u
            Q = sol(mids; idxs = sys.Q_in).u
            Qj = sol(mids; idxs = sys.jc.port_j.Q_flow).u
            R = Rh.(Ref(p), Tj)
            Qh = p.duty_on .* (p.V_sup ./ (p.R_load .+ R)) .^ 2 .* R
            @test all(Q[.!on] .== 0.0)                               # zero heat between bursts, exactly
            err = maximum(abs.(Q[on] ./ Qh[on] .- 1))
            worst = max(worst, err)
            @test err <= 1e-12
            @test maximum(abs.(Qj .- Q)) <= 1e-12 * maximum(Qh)      # Q_in is what enters the junction port
        end
        @info "ObservationHarness bursts: max rel |Q_in/(duty·I_d²·R_ds) − 1| in bursts = $(worst) (bar 1e-12); Q_in = 0 between"
    end

    @testset "case side: the measured s1 T_c, linearly interpolated" begin
        worst = 0.0
        for (nm, sol) in pairs(sols)
            Tc = sol(ts; idxs = sys.jc.port_c.T).u
            err = maximum(abs.(Tc .- Tc_hand.(ts)))
            worst = max(worst, err)
            @test err <= 1e-9
        end
        @info "ObservationHarness case side: max |T_c − s1 interpolant| = $(worst) K (bar 1e-9 K)"
    end

    @testset "junction energy balance: C_die·ΔT_j = ∫Q_in dt − ∫(T_j − T_c)/(R_jc(1 + a)) dt" begin
        # Integral identity over [0, stop]: 5-point Gauss–Legendre on every solver step (the steps include the burst
        # edges, so Q_in is smooth inside each), with Q_in and T_c computed by hand from T_j and the data.
        xg = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664]
        wg = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]
        worst = 0.0
        for (nm, sol) in pairs(sols)
            p = pars(sol)
            tk = sol.t
            h = diff(tk); c = (tk[1:end-1] .+ tk[2:end]) ./ 2
            tn = vec([c[k] + h[k] / 2 * xg[j] for j in 1:5, k in eachindex(h)])
            wn = vec([h[k] / 2 * wg[j] for j in 1:5, k in eachindex(h)])
            Tj = sol(tn; idxs = sys.T_j).u
            R = Rh.(Ref(p), Tj)
            Qin = ifelse.(burst.(tn), p.duty_on .* (p.V_sup ./ (p.R_load .+ R)) .^ 2 .* R, 0.0)
            Qout = (Tj .- Tc_hand.(tn)) ./ (p.R_jc * (1 + p.a))
            Ein = sum(wn .* Qin)
            lhs = p.C_die * (sol[sys.T_j][end] - sol[sys.T_j][1])
            rhs = Ein - sum(wn .* Qout)
            err = abs(lhs - rhs) / Ein
            worst = max(worst, err)
            @test Ein > 0
            @test err <= 1e-6
        end
        @info "ObservationHarness energy balance: max |C_die·ΔT_j − ∫(Q_in − Q)dt| / ∫Q_in dt = $(worst) (bar 1e-6)"
    end
end
