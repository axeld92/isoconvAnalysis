@testset "Recorded-time TG intervals and quadrature" begin
    mktempdir() do directory
        function recorded_curve(
            name,
            alpha,
            temperatures;
            times=collect(eachindex(alpha)),
            rows=eachindex(alpha),
        )
            source = joinpath(directory, "$name.csv")
            open(source, "w") do io
                println(io, times === nothing ? "T,m" : "T,t,m")
                for i in eachindex(alpha)
                    i == 3 && println(io)
                    mass = 8 * (1 - alpha[i])
                    if times === nothing
                        println(io, "$(temperatures[i]),$mass")
                    else
                        println(io, "$(temperatures[i]),$(times[i]),$mass")
                    end
                end
            end
            signal = ThermalSignalSpec(:mass, ThermalColumnSpec("m", :mg), :mass)
            mapping = ThermalColumnMapping(
                ThermalColumnSpec("T", :K),
                [signal];
                time=times === nothing ? nothing : ThermalColumnSpec("t", :min),
            )
            run = read_thermal_delimited(source, mapping; id=name, material_id="test")
            return tg_conversion(
                run,
                :mass;
                segment=1,
                reaction_rows=first(rows):last(rows),
                initial_reference=8,
                terminal_reference=0,
            )
        end

        alpha = [0, 0.25, 0.50, 0.75, 1]
        temperatures = [300, 400, 360, 490, 520]
        times = [0, 1, 3, 6, 10]
        curve = recorded_curve("recorded", alpha, temperatures; times)

        @testset "Clipped measurements and boundaries" begin
            interval = tg_recorded_time_interval(curve, 0.625; delta_alpha=0.5)
            @test interval.status === :supported
            @test interval.lower_crossing_status === :unique_upward
            @test interval.upper_crossing_status === :unique_upward
            @test interval.alpha == [0.125, 0.25, 0.5, 0.625]
            @test interval.time_min == [0.5, 1, 3, 4.5]
            @test interval.temperature_K == [350, 400, 360, 425]
            @test interval.source_lines == [(2, 3), (3, 3), (5, 5), (5, 6)]
            @test interval.segment_index == curve.segment_index
            @test interval.source_sha256 == curve.source_sha256
            @test interval.reversal_steps == 0
            @test interval.plateau_steps == 0
            @test interval.max_backtracking_alpha == 0
            @test interval.plateau_duration_min == 0

            exact = tg_recorded_time_interval(curve, 0.75; delta_alpha=0.5)
            @test exact.status === :supported
            @test exact.time_min == [1, 3, 6]
            @test exact.temperature_K == [400, 360, 490]
            @test exact.source_lines == [(3, 3), (5, 5), (6, 6)]

            same_step = tg_recorded_time_interval(curve, 0.4375; delta_alpha=0.0625)
            @test same_step.status === :supported
            @test same_step.time_min == [2, 2.5]
            @test same_step.temperature_K == [380, 370]
            @test same_step.source_lines == [(3, 5), (3, 5)]

            cropped = recorded_curve("cropped", alpha, temperatures; times, rows=2:4)
            supported = tg_recorded_time_interval(cropped, 0.625; delta_alpha=0.25)
            @test supported.status === :supported
            @test supported.source_lines == [(3, 5), (5, 5), (5, 6)]
            endpoint = tg_recorded_time_interval(cropped, 0.625; delta_alpha=0.375)
            @test endpoint.status === :boundary_not_unique
            @test endpoint.lower_crossing_status === :endpoint
            @test isempty(endpoint.time_min)
            missing = tg_recorded_time_interval(cropped, 0.625; delta_alpha=0.5)
            @test missing.lower_crossing_status === :missing
            for (upper, width) in
                ((0.05, 0.05), (0.01, 0.02), (1.0, 0.02), (0.5, eps() / 16))
                outside = tg_recorded_time_interval(curve, upper; delta_alpha=width)
                @test outside.status === :outside_interior
                @test outside.lower_crossing_status === :unchecked
            end
            @test_throws ArgumentError tg_recorded_time_interval(curve, 0.5; delta_alpha=0)
            @test_throws ArgumentError tg_recorded_time_interval(
                curve, 1.5; delta_alpha=0.02
            )
        end

        @testset "Whole-branch ambiguity and time support" begin
            reverse_alpha = [0, 0.125, 0.375, 0.50, 0.4375, 0.625, 0.875, 1]
            reverse_curve = recorded_curve("reversal", reverse_alpha, fill(450, 8))
            reverse = tg_recorded_time_interval(reverse_curve, 0.75; delta_alpha=0.5)
            @test reverse.lower_crossing_status === :unique_upward
            @test reverse.upper_crossing_status === :unique_upward
            @test reverse.status === :branch_reversal
            @test reverse.reversal_steps == 1
            @test reverse.max_backtracking_alpha == 0.0625
            @test_throws ArgumentError tg_recorded_time_quadrature(reverse, 60)

            flat_alpha = [0, 0.125, 0.375, 0.375, 0.625, 0.875, 1]
            flat_curve = recorded_curve("plateau", flat_alpha, fill(450, 7))
            flat = tg_recorded_time_interval(flat_curve, 0.75; delta_alpha=0.5)
            @test flat.status === :branch_plateau
            @test flat.plateau_steps == 1
            @test flat.plateau_duration_min == 1
            @test flat.max_backtracking_alpha == 0
            @test_throws ArgumentError tg_recorded_time_quadrature(flat, 60)

            ambiguous = tg_recorded_time_interval(reverse_curve, 0.75; delta_alpha=0.28125)
            @test ambiguous.status === :boundary_not_unique
            @test ambiguous.lower_crossing_status === :ambiguous
            @test ambiguous.upper_crossing_status === :unique_upward

            untimed = recorded_curve("untimed", alpha, temperatures; times=nothing)
            unsupported = tg_recorded_time_interval(untimed, 0.625; delta_alpha=0.5)
            @test unsupported.status === :missing_time
            @test unsupported.lower_crossing_status === :unique_upward
            @test isempty(unsupported.time_min)
            @test unsupported.plateau_steps === nothing
        end

        @testset "Fixed-energy quadrature checks" begin
            isothermal = recorded_curve("isothermal", alpha, fill(450, 5); times)
            interval = tg_recorded_time_interval(isothermal, 0.625; delta_alpha=0.5)
            isothermal_check = tg_recorded_time_quadrature(interval, 120)
            gas_constant = 0.00831446261815324
            expected = log(4.0) - 120 / (gas_constant * 450)
            for field in (
                :recorded_log_integral,
                :legacy_log_integral,
                :refined_log_integral,
                :doubled_log_integral,
            )
                @test getproperty(isothermal_check, field) ≈ expected atol=1e-12
            end

            thermal = tg_recorded_time_interval(curve, 0.625; delta_alpha=0.5)
            zero = tg_recorded_time_quadrature(thermal, 0)
            @test zero.recorded_log_integral ≈ log(4.0)
            @test zero.legacy_log_integral ≈ log(4.0)
            @test zero.doubled_log_integral ≈ log(4.0)
            probe = tg_recorded_time_quadrature(
                thermal, 60; legacy_alpha_nodes=3, subdivisions_per_step=32
            )
            @test probe.recorded_log_integral >
                probe.refined_log_integral >
                probe.doubled_log_integral
            @test abs(probe.legacy_log_integral - probe.doubled_log_integral) > 0.1
            t, T = thermal.time_min, thermal.temperature_K
            # Independent midpoint check on the declared piecewise-linear history.
            n = 4096
            midpoint = sum(2:length(t)) do j
                sum(1:n) do k
                    f = (k - 0.5) / n
                    temp = T[j - 1] + f * (T[j] - T[j - 1])
                    (t[j] - t[j - 1]) / n * exp(-60 / (gas_constant * temp))
                end
            end
            @test probe.doubled_log_integral ≈ log(midpoint) atol=0.001
            @test abs(probe.doubled_log_integral - log(midpoint)) <
                abs(probe.refined_log_integral - log(midpoint))
            @test_throws ArgumentError tg_recorded_time_quadrature(thermal, -1)
            @test_throws ArgumentError tg_recorded_time_quadrature(
                thermal, 60; legacy_alpha_nodes=2
            )
            @test_throws ArgumentError tg_recorded_time_quadrature(
                thermal, 60; subdivisions_per_step=1
            )
        end
    end
end
