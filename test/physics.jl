# Physics checks of the Stage A Dyad components (dyad/*.dyad) against hand calculations.
using PowerElectronics, DyadInterface, ModelingToolkit, Test

R_on(m, T) = m.R₀ * (1 + m.α₁ * (T - m.T_ref) + m.α₂ * (T - m.T_ref)^2)           # pristine (ΔR = 0)
const MOS = (R₀ = 0.13210477776821447, T_ref = 298.15, α₁ = 0.004971635253450363, α₂ = 3.121126781011276e-5)   # ConductionMOSFET defaults
const RIG = (V_sup = 5.61, R_load = 0.2, duty = 0.4, T_amb = 298.15, R_ca = 58.9, C_case = 5.0,
             R_th = [6.634289458535688, 4.058275425543876], C_th = [0.03132956205116243, 0.6933552281698974])   # NASARig / CauerNetwork defaults

@testset "NASARig switching window" begin
    res = NASARigSwitching()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    # mid on-pulse (period 1 ms, on for 0.4 ms): KVL around supply -> load -> channel
    for t in (0.2e-3, 2.2e-3, 4.2e-3)
        Tj = sol(t, idxs = m.T_j)
        @test sol(t, idxs = m.mosfet.i) ≈ RIG.V_sup / (RIG.R_load + R_on(MOS, Tj)) rtol = 1e-6
    end
    # mid off-interval: leakage only, the drain sits at V_sup
    @test abs(sol(0.7e-3, idxs = m.mosfet.i)) < 1e-4
    @test sol(0.7e-3, idxs = m.mosfet.v) ≈ RIG.V_sup rtol = 1e-4
    # heat delivered over 5 periods ≈ duty × on-state power (T_j changes by a few K within the window)
    ts = range(0, 5e-3, length = 50_001)
    P = sol(ts, idxs = m.mosfet.P_loss).u
    E_sw = sum(P) * step(ts)
    Tj0 = RIG.T_amb; i_on = RIG.V_sup / (RIG.R_load + R_on(MOS, Tj0))
    E_avg = RIG.duty * R_on(MOS, Tj0) * i_on^2 * 5e-3
    @test E_sw ≈ E_avg rtol = 0.03
end

@testset "AveragedRig warm-up: energy balance" begin
    res = AveragedRigWarmup()
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    # on-state operating point at every saved time: i_on = V_sup / (R_load + R_on(T_j))
    Tj = sol[m.T_j]; i = sol[m.mosfet.i]
    @test maximum(abs.(i .- RIG.V_sup ./ (RIG.R_load .+ R_on.(Ref(MOS), Tj)))) < 1e-6
    # energy in = stored + lost to ambient
    ts = range(0, 600, length = 60_001); dt = step(ts)
    Pin = sol(ts, idxs = m.P_avg).u; Tc = sol(ts, idxs = m.T_c).u
    E_in = sum(Pin) * dt
    E_out = sum((Tc .- RIG.T_amb) ./ RIG.R_ca) * dt
    T1, T2 = sol(600.0, idxs = m.jc.T[1]), sol(600.0, idxs = m.jc.T[2])
    E_store = RIG.C_th[1] * (T1 - RIG.T_amb) + RIG.C_th[2] * (T2 - RIG.T_amb) + RIG.C_case * (Tc[end] - RIG.T_amb)
    @test E_in ≈ E_store + E_out rtol = 1e-3
end

@testset "AveragedRig steady state" begin
    duty = 0.05
    model = PowerElectronics.AveragedRig(; name = :AveragedRig, duty = duty)
    res = AveragedRigWarmup(; model = model, stop = 5000.0)        # 17 case time constants (R_ca·C_case ≈ 295 s)
    sol, m = res.sol, symbolic_container(res)
    @test string(sol.retcode) == "Success"
    # fixed point T_j = T_amb + P(T_j) (ΣR_th + R_ca), P = duty · R_on · (V_sup / (R_load + R_on))²
    Rtot = sum(RIG.R_th) + RIG.R_ca
    Tj = RIG.T_amb
    for _ in 1:200
        R = R_on(MOS, Tj); Tj = RIG.T_amb + duty * R * (RIG.V_sup / (RIG.R_load + R))^2 * Rtot
    end
    @test sol[m.T_j][end] ≈ Tj rtol = 1e-4
    @test sol[m.T_c][end] ≈ RIG.T_amb + sol[m.P_avg][end] * RIG.R_ca rtol = 1e-4
end
