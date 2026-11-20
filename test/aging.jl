using PowerElectronics, DyadInterface, ModelingToolkit, Test

@testset "Aging input and AgingState" begin
    @testset "JunctionStage: a_in held at 0.5 gives the solution of the parameter a = 0.5, bit for bit" begin
        par = JunctionStep(model = PowerElectronics.JunctionStepTest(; name = :JunctionStepTest, a_net = 0.5))
        sig = JunctionStepInput(model = PowerElectronics.JunctionStepInputTest(; name = :JunctionStepInputTest, a_net = 0.5))
        mp, ms = symbolic_container(par), symbolic_container(sig)
        @test sig.sol.t == par.sol.t
        @test sig.sol[ms.jc.T_j] == par.sol[mp.jc.T_j]
        @test sig.sol[ms.jc.T_j][end] - 400.0 ≈ 2.0 * 3.1 * 1.5 rtol = 1e-8
    end

    @testset "CauerNetwork: a_in held at 0.5 gives the solution of the parameter a = 0.5, bit for bit" begin
        par = CauerStep(model = PowerElectronics.CauerStepTest(; name = :CauerStepTest, a_net = 0.5))
        sig = CauerStepInput(model = PowerElectronics.CauerStepInputTest(; name = :CauerStepInputTest, a_net = 0.5))
        mp, ms = symbolic_container(par), symbolic_container(sig)
        @test sig.sol.t == par.sol.t
        @test sig.sol[ms.jc.T[1]] == par.sol[mp.jc.T[1]]
        @test sig.sol[ms.jc.T[2]] == par.sol[mp.jc.T[2]]
        # aging on stage 1 only: P·(R_th[1]·(1 + a) + R_th[2]) = 2·(1.0·1.5 + 2.1) = 7.2 K at node 1 and
        # P·R_th[2] = 2·2.1 = 4.2 K at node 2 (both stages scaled would give 9.3 K and 6.3 K)
        @test sig.sol[ms.jc.T[1]][end] - 300.0 ≈ 2.0 * (1.0 * 1.5 + 2.1) rtol = 1e-8
        @test sig.sol[ms.jc.T[2]][end] - 300.0 ≈ 2.0 * 2.1 rtol = 1e-8
    end

    @testset "CauerNetwork with N = 1: the single stage is scaled, and a_in held at 0.5 gives the solution of the parameter a = 0.5, bit for bit" begin
        one_stage = Dict{String, Any}("jc.N" => 1, "jc.R_th" => [1.0], "jc.C_th" => [0.005])
        par = CauerStep(model = PowerElectronics.CauerStepTest(; name = :CauerStepTest, a_net = 0.5, __overrides = one_stage))
        sig = CauerStepInput(model = PowerElectronics.CauerStepInputTest(; name = :CauerStepInputTest, a_net = 0.5, __overrides = one_stage))
        mp, ms = symbolic_container(par), symbolic_container(sig)
        @test sig.sol.t == par.sol.t
        @test sig.sol[ms.jc.T[1]] == par.sol[mp.jc.T[1]]
        # one stage of 1.0 K/W: P·R_th[1]·(1 + a) = 2·1.0·1.5 = 3.0 K
        @test sig.sol[ms.jc.T[1]][end] - 300.0 ≈ 2.0 * 1.0 * 1.5 rtol = 1e-8
    end

    @testset "ConductionMOSFET: ΔR_in held at 0.05 Ω gives the solution of the parameter ΔR = 0.05 Ω, bit for bit" begin
        par = MOSFETConduction()
        sig = MOSFETConductionInput()
        mp, ms = symbolic_container(par), symbolic_container(sig)
        @test sig.sol[ms.mosfet.i] == par.sol[mp.mosfet.i]
        @test sig.sol[ms.mosfet.R_on] == par.sol[mp.mosfet.R_on]
        @test sig.sol[ms.mosfet.R_on][end] ≈ 0.18 * (1 + 0.0117 * 100 + 2e-5 * 100^2) + 0.05 rtol = 1e-12
    end

    @testset "AgingState: a law follows its closed form and acts on R_jc while it grows" begin
        res = AgingClosedForm()
        sol, m = res.sol, symbolic_container(res)
        P_in, T_case, R_jc = sol.ps[m.P_in], sol.ps[m.T_case], sol.ps[m.jc.R_jc]
        a_start, a_s = sol.ps[m.a_start], sol.ps[m.law.a_s]
        κ_T = sol.ps[m.law.κ] * exp(-sol.ps[m.law.θ] * (1 / T_case - 1 / sol.ps[m.law.T_law]))
        a = sol[m.law.a]
        @test a[1] == a_start
        @test a[end] > 7 * a_start
        @test maximum(abs.(a ./ ((a_start + a_s) .* exp.(κ_T .* sol.t) .- a_s) .- 1)) < 1e-8
        @test maximum(abs.(sol[m.jc.Q] .* R_jc .* (1 .+ a) .- (sol[m.jc.T_j] .- T_case))) < 1e-9
        @test sol[m.jc.T_j][end] - T_case ≈ P_in * R_jc * (1 + a[end]) rtol = 1e-4
    end
end
