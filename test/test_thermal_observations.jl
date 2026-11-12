@testset "thermal observations and segment selection" begin
    mass = ThermalSignalSpec(:tg, ThermalColumnSpec("mass", :mg), :mass)
    dsc = ThermalSignalSpec(
        :dsc, ThermalColumnSpec("heat", :mW_per_mg), :specific_heat_flow; exotherm_sign=-1
    )
    aux = ThermalSignalSpec(:voltage, ThermalColumnSpec("voltage", :uV), :auxiliary)
    mapping = ThermalColumnMapping(
        ThermalColumnSpec("temp", :C),
        (mass, dsc, aux);
        time=ThermalColumnSpec("time", :s),
        segment_column="phase",
    )
    mktempdir() do directory
        path = joinpath(directory, "measured.csv")
        header = "temp,time,mass,heat,voltage,phase\n"
        data =
            header *
            "30,60,20,-1,4,heat\n40,120,19,-2,5,heat\n" *
            "40,0,18,0,6,hold\n40,60,17,1,7,hold\n" *
            "40,0,16,2,8,heat\n35,60,15,3,9,heat\n30,0,14,4,10,end\n"
        write(path, data)
        run = read_thermal_delimited(path, mapping)
        @test temperature_K(run) ≈ [303.15, 313.15, 313.15, 313.15, 313.15, 308.15, 303.15]
        @test time_min(run) == [1, 2, 0, 1, 0, 1, 0]
        @test signal_values(run, :tg; unit=:g) ≈
            [0.020, 0.019, 0.018, 0.017, 0.016, 0.015, 0.014]
        @test signal_values(run, :dsc; unit=:W_per_kg) ==
            [-1000, -2000, 0, 1000, 2000, 3000, 4000]
        @test signal_values(run, :voltage; unit=:uV) == collect(4:10)
        @test_throws ArgumentError signal_values(run, :tg; unit=:fraction)
        @test_throws ArgumentError signal_values(run, :dsc; unit=:W)
        @test_throws ArgumentError signal_values(run, :voltage; unit=:mW)
        @test_throws ArgumentError signal_values(run, :absent)

        diagnostics = thermal_segment_diagnostics(run)
        @test getproperty.(diagnostics, :slope_K_per_min) == [10, 0, -5, nothing]
        @test getproperty.(diagnostics, :duration_min) == [1, 1, 1, 0]
        @test getproperty.(diagnostics, :increasing_steps) == [1, 0, 0, 0]
        @test getproperty.(diagnostics, :decreasing_steps) == [0, 0, 1, 0]
        @test getproperty.(diagnostics, :flat_steps) == [0, 1, 0, 0]
        selection = select_thermal_segment(run, 3)
        segment = only(selection.segments)
        @test (segment.index, segment.label, segment.rows) == (3, "heat", 1:2)
        @test selection.provenance.source_lines == [6, 7]
        @test selection.provenance.sha256 == run.provenance.sha256
        @test selection.time == [0, 60]
        @test selection.signals[2].spec.exotherm_sign == -1
        @test signal_values(selection, :dsc) == [2, 3]
        @test_throws ArgumentError select_thermal_segment(selection, 1)
        selection.signals[2].values[1] = 99
        @test signal_values(run, :dsc)[5] == 2

        write(path, replace(data, "40,120,19" => "40,60,19"))
        @test_throws r"time must strictly increase within segment 1 at source line 3" read_thermal_delimited(
            path, mapping
        )
        write(path, replace(data, "30,60,20" => "-273.15,60,20"))
        @test_throws r"above absolute zero at source line 2" read_thermal_delimited(
            path, mapping
        )

        relative = ThermalSignalSpec(:tg, ThermalColumnSpec("mass", :percent), :mass)
        no_time = ThermalColumnMapping(ThermalColumnSpec("temp", :K), [relative])
        write(path, "temp,mass\n300,99\n350,93\n")
        unknown_program = read_thermal_delimited(path, no_time)
        @test only(unknown_program.segments).label === nothing
        @test temperature_K(unknown_program) == [300, 350]
        @test time_min(unknown_program) === nothing
        @test signal_values(unknown_program, :tg; unit=:fraction) ≈ [0.99, 0.93]
        diagnostic = only(thermal_segment_diagnostics(unknown_program))
        @test diagnostic.duration_min === nothing
        @test diagnostic.slope_K_per_min === nothing
        @test diagnostic.increasing_steps == 1
    end

    @test_throws ArgumentError ThermalColumnMapping(ThermalColumnSpec("temp", :F), [mass])
    @test_throws ArgumentError ThermalColumnMapping(ThermalColumnSpec("mass", :C), [mass])
    @test_throws ArgumentError ThermalColumnMapping(
        ThermalColumnSpec("temp", :C), [mass, mass]
    )
    @test_throws ArgumentError ThermalSignalSpec(:tg, mass.column, :mass; exotherm_sign=-1)
    @test_throws ArgumentError ThermalSignalSpec(
        :dsc, dsc.column, :specific_heat_flow; exotherm_sign=0
    )
    @test_throws ArgumentError ThermalSignalSpec(:dsc, dsc.column, :heat_flow)
end
