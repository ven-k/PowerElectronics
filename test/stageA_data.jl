# Stage A input files (assets/stageA/, derived from the NASA data): every file listed in
# index.csv loads as a DyadTimeseries through its dyad:// URI, as a Dyad analysis reads it, and keeps the
# invariants the calibration relies on.
using DyadData, CSV, Test

const STAGEA = joinpath(pkgdir(PowerElectronics), "assets", "stageA")

@testset "Stage A input files" begin
    index = CSV.File(joinpath(STAGEA, "index.csv"))
    @test length(index) > 0
    for r in index
        if r.kind == "bursts"                                         # event-resolution heat input
            tbl = DyadTable("dyad://PowerElectronics/stageA/$(r.file)"; columns = ["on", "off", "P"])
            on, off, P = tbl["on"], tbl["off"], tbl["P"]
            @testset "$(r.file)" begin
                @test length(on) == r.rows
                @test all(off .> on) && all(on[2:end] .>= off[1:end-1])   # ordered, no overlap
                @test all(isfinite, P) && all(P .>= 0)
            end
            continue
        end
        if r.kind == "heat"                                           # step series for ConstantInterpolation
            ts = DyadTimeseries("dyad://PowerElectronics/stageA/$(r.file)"; independent_var = "time", dependent_vars = ["P"])
            t, P = ts["time"], ts["P"]
            b = CSV.File(joinpath(STAGEA, "test$(r.test)_bursts.csv"))
            E_burst = sum(row.P * (row.off - row.on) for row in b)
            E_step = sum(P[i] * (t[i+1] - t[i]) for i in 1:length(t)-1)
            @testset "$(r.file)" begin
                @test length(t) == r.rows
                @test t[1] == 0.0                                     # a burst at time 0 carries its P in row 1
                @test all(diff(t) .> 0)                               # coincident edges deduplicated
                @test !any(ismissing, P) && all(isfinite, P) && all(P .>= 0)
                @test abs(E_step - E_burst) <= 1e-9 * E_burst         # the square wave carries the burst energy
            end
            continue
        end

        cols = r.kind == "resistance" ? ["R_meas", "T_c", "phase", "baseline", "t_burst"] : ["P_heat", "T_c", "switching", "phase"]
        ts = DyadTimeseries("dyad://PowerElectronics/stageA/$(r.file)"; independent_var = "time", dependent_vars = cols)
        t = ts["time"]
        @testset "$(r.file)" begin
            @test length(t) == r.rows
            @test t[1] >= 0
            @test all(diff(t) .> 0)                                   # strictly increasing: a valid ODE input grid
            for c in cols
                @test !any(ismissing, ts[c]) && all(isfinite, ts[c])  # no gaps: CSV select would load them as missing
            end
            @test all(253.15 .<= ts["T_c"] .<= 673.15)                # -20 ... 400 C, kelvin
            @test all(in((1, 2, 3)), ts["phase"])
            if r.kind == "resistance"
                @test all(0.05 .<= ts["R_meas"] .<= 3.0)
                @test all(in((0, 1)), ts["baseline"]) && sum(ts["baseline"]) > 0
                @test all(ts["t_burst"] .>= -1)
            else
                @test all(ts["P_heat"] .>= 0)
                @test all(0 .<= ts["switching"] .<= 1)
                r.kind == "cooldown" && @test all(==(0), ts["P_heat"]) && all(==(3), ts["phase"])
            end
        end
    end
end
