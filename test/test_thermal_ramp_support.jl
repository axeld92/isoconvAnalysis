@testset "Recorded TG constant-ramp diagnostics" begin
    mktempdir() do directory
        function ramp_curve(
            name, temperature, progress; clock=nothing, rows=eachindex(progress)
        )
            path = joinpath(directory, "$name.csv")
            open(path, "w") do file
                println(
                    file, clock === nothing ? "temperature,mass" : "temperature,mass,time"
                )
                for i in eachindex(progress)
                    if clock === nothing
                        println(file, temperature[i], ',', 8 * (1 - progress[i]))
                    else
                        println(
                            file, temperature[i], ',', 8 * (1 - progress[i]), ',', clock[i]
                        )
                    end
                end
            end
            mapping = ThermalColumnMapping(
                ThermalColumnSpec("temperature", :K),
                [ThermalSignalSpec(:mass, ThermalColumnSpec("mass", :mg), :mass)];
                time=clock === nothing ? nothing : ThermalColumnSpec("time", :min),
            )
            run = read_thermal_delimited(path, mapping; material_id="ramp control")
            return tg_conversion(
                run,
                :mass;
                segment=1,
                reaction_rows=first(rows):last(rows),
                initial_reference=8,
                terminal_reference=0,
            )
        end

        linear = ramp_curve(
            "linear", [300, 340, 380, 400], [0, 0.25, 0.75, 1]; clock=[0, 2, 4, 5]
        )
        result = tg_ramp_diagnostic(
            linear, 0.50; start_row=1, probe_energies_kJ_per_mol=[60]
        )
        @test result.status === :evaluated
        @test result.start_source_line == 2
        @test result.raw_row_count == 2
        @test result.duration_min ≈ 3.0
        @test result.secant_K_per_min ≈ 20
        @test result.OLS_K_per_min ≈ 20
        @test collect(result.half_secants_K_per_min) ≈ [20, 20]
        @test result.affine_max_abs_residual_K < 1e-12
        @test result.falling_temperature_steps == 0
        probe = only(result.probes)
        @test probe.measured_vs_affine_log_ratio ≈ 0 atol=1e-13
        # Independent adaptive quadrature of exp(-60/(R*(300+20t))) over 0–3 min.
        @test exp(probe.log_measured_integral_min) ≈ 1.5937538825291003e-9 rtol=1e-12
        @test probe.max_quadrature_log_difference < 1e-5
        @test 0 < probe.omitted_start_fraction < 1

        nonlinear = ramp_curve(
            "nonlinear",
            [300, 350, 340, 380, 400],
            [0, 0.125, 0.25, 0.75, 1];
            clock=[0, 1, 2, 4, 5],
        )
        departure = tg_ramp_diagnostic(
            nonlinear, 0.50; start_row=1, probe_energies_kJ_per_mol=[60, 120]
        )
        @test departure.status === :evaluated
        @test departure.secant_K_per_min ≈ 20
        @test departure.affine_max_abs_residual_K ≈ 30
        @test departure.falling_temperature_steps == 1
        @test departure.cumulative_cooling_K ≈ 10
        @test collect(departure.step_slope_range_K_per_min) ≈ [-10, 50]
        @test all(p -> p.measured_vs_affine_log_ratio > 0, departure.probes)
        warm = tg_ramp_diagnostic(
            nonlinear, 0.50; start_row=3, probe_energies_kJ_per_mol=[60]
        )
        @test warm.duration_min ≈ 1.0
        @test warm.start_source_line == 4
        @test only(warm.probes).omitted_start_fraction > probe.omitted_start_fraction
        @test only(warm.probes).measured_vs_affine_log_ratio ≈ 0 atol=1e-13

        crop = ramp_curve(
            "crop", [300, 340, 380, 400], [0, 0.25, 0.75, 1]; clock=[0, 2, 4, 5], rows=2:4
        )
        @test tg_ramp_diagnostic(crop, 0.50; start_row=2).start_time_min == 2
        @test tg_ramp_diagnostic(crop, 0.50; start_row=3).status ===
            :start_at_or_after_crossing
        @test tg_ramp_diagnostic(crop, 0.10; start_row=2).status === :unsupported_crossing
        @test_throws ArgumentError tg_ramp_diagnostic(crop, 0.50; start_row=1)
        @test_throws ArgumentError tg_ramp_diagnostic(
            linear, 0.5; start_row=1, probe_energies_kJ_per_mol=[0]
        )
        @test_throws ArgumentError tg_ramp_diagnostic(
            linear, 0.5; start_row=1, probe_energies_kJ_per_mol=[60, 60]
        )

        untimed = ramp_curve("untimed", [300, 340, 380], [0, 0.50, 1])
        @test tg_ramp_diagnostic(untimed, 0.50; start_row=1).status === :missing_time
        cooling = ramp_curve("cooling", [360, 340, 320], [0, 0.50, 1]; clock=[0, 1, 2])
        reverse = tg_ramp_diagnostic(
            cooling, 0.50; start_row=1, probe_energies_kJ_per_mol=[60]
        )
        @test reverse.status === :nonpositive_net_ramp
        @test reverse.falling_temperature_steps == 1
        @test reverse.duration_min == 1
        @test isempty(reverse.probes)
    end

    numerical = IsoconversionalAnalysis._ramp_log_temperature_integral
    @test exp(numerical(20.0)) ≈ 4.7024282154283385e-12 rtol=1e-10
    @test exp(numerical(0.10)) ≈ 7.225450221940205 rtol=1e-10
    @test exp(numerical(1.0)) ≈ 0.14849550677592205 rtol=1e-10
end
