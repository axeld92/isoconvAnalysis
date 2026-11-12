@testset "TG local Friedman rate evidence" begin
    mktempdir() do directory
        function rate_curve(name, times, alpha; recorded_time=true)
            path = joinpath(directory, "$name.csv")
            open(path, "w") do io
                println(io, recorded_time ? "temperature,time,mass" : "temperature,mass")
                for (t, a) in zip(times, alpha)
                    temperature = 300 + 5t
                    mass = 8 * (1 - a)
                    println(
                        io, recorded_time ? "$temperature,$t,$mass" : "$temperature,$mass"
                    )
                end
            end
            mapping = ThermalColumnMapping(
                ThermalColumnSpec("temperature", :K),
                [ThermalSignalSpec(:m, ThermalColumnSpec("mass", :mg), :mass)];
                time=recorded_time ? ThermalColumnSpec("time", :min) : nothing,
            )
            run = read_thermal_delimited(path, mapping; material_id="polynomial_fixture")
            return tg_conversion(
                run,
                :m;
                segment=1,
                reaction_rows=1:length(alpha),
                initial_reference=8,
                terminal_reference=0,
            )
        end

        time = [0.0, 0.43, 1.2, 1.81, 2.3, 2.79, 3.4, 4.8]
        alpha = @. 0.5 + 0.1 * (time - 2) + 0.02 * (time - 2)^2
        curve = rate_curve("quadratic", time, alpha)
        sample = tg_friedman_rate(curve, 0.5; half_window_min=1.1, degree=2)
        @test sample.status === :positive_rate
        @test sample.crossing_status === :unique_upward
        @test sample.crossing_source_lines == (5, 6)
        @test sample.crossing_temperature_K ≈ 300 + 5sample.crossing_time_min
        @test sample.window_source_lines == (4, 7)
        @test sample.window_time_min == (1.2, 2.79)
        @test all(isapprox.(sample.window_conversion_range, extrema(alpha[3:6])))
        @test (sample.points_before, sample.points_after, sample.window_points) == (2, 2, 4)
        @test sample.conversion_rate_min_inv ≈ 0.1 + 0.04 * (sample.crossing_time_min - 2)
        @test sample.residual_rmse < 1e-14
        @test curve.alpha ≈ alpha
        @test curve.time_min == time

        narrow = tg_friedman_rate(curve, 0.5; half_window_min=0.4)
        @test narrow.status === :insufficient_window
        @test narrow.window_source_lines == (5, 6)
        @test narrow.window_points == 2
        @test narrow.conversion_rate_min_inv === nothing
        @test tg_friedman_rate(curve, 0.5; half_window_min=1.1, degree=3).status ===
            :insufficient_window
        @test tg_friedman_rate(curve, 0.5; half_window_min=1e300).status ===
            :fit_unavailable
        @test tg_friedman_rate(curve, 0.02; half_window_min=1.0).crossing_status ===
            :missing

        linear_time = collect(0.0:5.0)
        linear_alpha = @. 0.25 + 0.125 * linear_time
        linear = rate_curve("linear", linear_time, linear_alpha)
        exact = tg_friedman_rate(linear, 0.5; half_window_min=2.0, degree=3)
        @test exact.crossing_time_min == 2
        @test exact.window_points == 5
        @test exact.conversion_rate_min_inv ≈ 0.125
        @test tg_friedman_rate(linear, 0.25; half_window_min=2.0).status ===
            :crossing_unavailable
        untimed = rate_curve("untimed", linear_time, linear_alpha; recorded_time=false)
        @test tg_friedman_rate(untimed, 0.5; half_window_min=2.0).status ===
            :time_unavailable

        oscillatory = rate_curve(
            "oscillatory",
            [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            [0.0, 0.47, 0.46, 0.52, 1.0, 0.51],
        )
        negative = tg_friedman_rate(oscillatory, 0.5; half_window_min=2, degree=3)
        @test negative.crossing_status === :unique_upward
        @test negative.status === :nonpositive_rate
        @test negative.conversion_rate_min_inv ≈ -0.2390118064031843
        @test tg_friedman_rate(oscillatory, 0.5; half_window_min=2, degree=1).status ===
            :positive_rate
        @test_throws ArgumentError tg_friedman_rate(curve, 0.5; half_window_min=0)
        @test_throws ArgumentError tg_friedman_rate(curve, 0.5; half_window_min=Inf)
        @test_throws ArgumentError tg_friedman_rate(curve, 0.5; half_window_min=1, degree=4)
        @test_throws ArgumentError tg_friedman_rate(
            curve, 0.5; half_window_min=1, minimum_points_per_side=0
        )
    end
end
