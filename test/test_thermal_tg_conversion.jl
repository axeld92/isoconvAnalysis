@testset "TG conversion from measured mass" begin
    mass = ThermalSignalSpec(:m, ThermalColumnSpec("Mass", :mg), :mass)
    mapping = ThermalColumnMapping(
        ThermalColumnSpec("T", :C),
        [mass];
        time=ThermalColumnSpec("Clock", :s),
        segment_column="Phase",
    )
    mktempdir() do directory
        source = joinpath(directory, "remaining_mass.csv")
        write(
            source,
            """T,Clock,Mass,Phase
30,100,18,warm
50,110,17,warm
310,0,10,heat
330,60,9,heat

325,120,8,heat
380,180,8.5,heat
430,240,4,heat
490,300,1,heat
500,360,0.9,heat
""",
        )
        run = read_thermal_delimited(source, mapping; material_id="test material")
        result = tg_conversion(
            run, :m; segment=2, reaction_rows=2:6, initial_reference=8, terminal_reference=2
        )
        @test result.run_id == "remaining_mass"
        @test result.material_id == "test material"
        @test result.segment_index == 2
        @test result.segment_label == "heat"
        @test result.signal_id == :m
        @test result.mass_unit == :mg
        @test result.reaction_rows == 2:6
        @test result.source_path == abspath(source)
        @test result.source_sha256 == run.provenance.sha256
        @test result.source_lines == [5, 7, 8, 9, 10]
        @test result.temperature_K ≈ [330, 325, 380, 430, 490] .+ 273.15
        @test result.time_min == [1, 2, 3, 4, 5]
        @test result.mass == [9, 8, 8.5, 4, 1]
        @test result.alpha ≈ [-1 / 6, 0, -1 / 12, 2 / 3, 7 / 6]
        @test result.decreasing_steps == 1
        @test result.outside_unit_interval_count == 3
        @test result.negative_mass_count == 0
        @test result.initial_reference.method == :fixed
        @test result.initial_reference.segment_rows === nothing
        @test isempty(result.initial_reference.source_lines)
        @test result.initial_reference.observed_mass_range === nothing

        medians = tg_conversion(
            select_thermal_segment(run, 2),
            :m;
            segment=2,
            reaction_rows=2:6,
            initial_reference=1:2,
            terminal_reference=6:7,
        )
        @test medians.source_lines == result.source_lines
        @test medians.mass == result.mass
        @test medians.initial_reference.mass == 9.5
        @test medians.initial_reference.method == :median
        @test medians.initial_reference.segment_rows == 1:2
        @test medians.initial_reference.source_lines == [4, 5]
        @test medians.initial_reference.observed_mass_range == (9.0, 10.0)
        @test medians.terminal_reference.mass ≈ 0.95
        @test medians.terminal_reference.source_lines == [10, 11]
        @test medians.alpha ≈ (9.5 .- result.mass) ./ (9.5 - 0.95)
        medians.mass[1] = 100
        @test run.signals[1].values[4] == 9

        options = (
            segment=2, reaction_rows=2:6, initial_reference=1:2, terminal_reference=6:7
        )
        @test_throws ArgumentError tg_conversion(run, :missing; options...)
        @test_throws ArgumentError tg_conversion(run, :m; options..., segment=3)
        @test_throws ArgumentError tg_conversion(run, :m; options..., reaction_rows=0:2)
        @test_throws ArgumentError tg_conversion(run, :m; options..., reaction_rows=4:3)
        @test_throws ArgumentError tg_conversion(
            run, :m; options..., terminal_reference=6:8
        )
        @test_throws ArgumentError tg_conversion(run, :m; options..., initial_reference=1:6)
        @test_throws ArgumentError tg_conversion(run, :m; options..., initial_reference=NaN)
        @test_throws ArgumentError tg_conversion(
            run, :m; options..., terminal_reference=Inf
        )
        @test_throws ArgumentError tg_conversion(run, :m; options..., terminal_reference=10)
    end

    mktempdir() do directory
        source = joinpath(directory, "apparent_mass.csv")
        write(
            source,
            "T,Mass,OtherMass,Loss\n300,10,1.1,0\n320,5,0.6,-5\n450,-0.2,0.08,-10.2\n",
        )
        grams = ThermalSignalSpec(:g, ThermalColumnSpec("OtherMass", :g), :mass)
        loss = ThermalSignalSpec(:loss, ThermalColumnSpec("Loss", :mg), :mass_loss)
        no_time = ThermalColumnMapping(ThermalColumnSpec("T", :K), [mass, grams, loss])
        run = read_thermal_delimited(source, no_time)
        options = (
            segment=1, reaction_rows=1:3, initial_reference=1:1, terminal_reference=3:3
        )
        apparent = tg_conversion(run, :m; options...)
        @test apparent.time_min === nothing
        @test apparent.terminal_reference.mass == -0.2
        @test apparent.negative_mass_count == 1
        @test apparent.alpha ≈ [0, 5 / 10.2, 1]
        @test apparent.alpha ≈ tg_conversion(run, :g; options...).alpha
        zero_reference = tg_conversion(run, :m; options..., terminal_reference=0)
        @test zero_reference.alpha ≈ [0, 0.5, 1.02]
        @test zero_reference.outside_unit_interval_count == 1
        @test_throws ArgumentError tg_conversion(run, :loss; options...)
    end
end
