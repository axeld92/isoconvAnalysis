@testset "Apparent energy from observed TG curves" begin
    mktempdir() do directory
        known_energy = 137.0
        gas = 0.00831446261815324
        temperatures = [495.0, 545.0, 610.0]
        exact_rates = 1e12 .* exp.(-known_energy ./ (gas .* temperatures))
        mapping = ThermalColumnMapping(
            ThermalColumnSpec("temperature", :K),
            [ThermalSignalSpec(:mass, ThermalColumnSpec("mass", :mg), :mass)];
            time=ThermalColumnSpec("time", :min),
        )
        function fixture(id, temperature, rate; maximum_alpha=1.0)
            path = joinpath(directory, "$id.csv")
            alpha = range(0.0, maximum_alpha; length=101)
            open(path, "w") do io
                println(io, "temperature,time,mass")
                for a in alpha
                    println(io, "$temperature,$(a / rate),$(9 - 8a)")
                end
            end
            run = read_thermal_delimited(path, mapping; id, material_id="analytic_tg")
            return tg_conversion(
                run,
                :mass;
                segment=1,
                reaction_rows=1:length(alpha),
                initial_reference=9,
                terminal_reference=1,
            )
        end
        curves = [
            fixture("run_$i", T, rate) for
            (i, (T, rate)) in enumerate(zip(temperatures, exact_rates))
        ]
        targets = [0.1, 0.5, 0.9]
        half_windows = 0.08 ./ exact_rates
        friedman = tg_activation_energy(
            :friedman, curves, targets; friedman_half_windows_min=half_windows
        )
        @test friedman.run_ids == ["run_1", "run_2", "run_3"]
        @test friedman.material_id == "analytic_tg"
        for point in friedman.points
            @test point.status === :estimated
            @test point.energy_kJ_per_mol ≈ known_energy atol=1e-9
            @test point.method_inputs ≈ exact_rates
            @test point.regression.r_squared ≈ 1.0
        end

        settings(delta, bounds=(20.0, 300.0)) =
            AdvancedVyazovkinConfig(delta, 41, bounds, 61, 1e-7, 200, 0.5, 1e-8)
        recorded = tg_activation_energy(
            :advanced_vyazovkin, curves, targets; advanced=settings(0.02)
        )
        for point in recorded.points
            @test point.status === :estimated
            @test point.energy_kJ_per_mol ≈ known_energy atol=1e-5
            @test point.objective ≈ 6.0 atol=1e-12
            @test point.run_statuses == fill(:supported, 3)
        end
        limited = only(
            tg_activation_energy(
                :advanced_vyazovkin, curves, [0.5]; advanced=settings(0.02, (40.0, 100.0))
            ).points,
        )
        @test limited.status === :boundary_solution
        @test limited.energy_kJ_per_mol ≈ 100.0 atol=1e-5
        unsupported = only(
            tg_activation_energy(
                :advanced_vyazovkin, curves, [0.1]; advanced=settings(0.2)
            ).points,
        )
        @test unsupported.status === :unsupported_interval
        @test unsupported.energy_kJ_per_mol === nothing

        for (method, power, coefficient) in
            ((:kas, 2.0, 1.0), (:fwo, 0.0, 1.052), (:starink, 1.92, 1.0008))
            beta = @. temperatures^power *
                exp(-coefficient * known_energy / (gas * temperatures))
            fit = only(
                tg_activation_energy(method, curves, [0.5]; ramp_rates_K_per_min=beta).points,
            )
            @test fit.status === :estimated
            @test fit.energy_kJ_per_mol ≈ known_energy atol=1e-9
        end

        narrow = only(
            tg_activation_energy(
                :friedman, curves, [0.5]; friedman_half_windows_min=half_windows ./ 100
            ).points,
        )
        @test narrow.status === :unsupported_rate
        @test narrow.energy_kJ_per_mol === nothing
        incomplete = [
            curves[1:2];
            fixture("incomplete", temperatures[3], exact_rates[3]; maximum_alpha=0.4)
        ]
        blocked = only(
            tg_activation_energy(
                :friedman, incomplete, [0.5]; friedman_half_windows_min=half_windows
            ).points,
        )
        @test blocked.status === :unsupported_crossing
        @test blocked.energy_kJ_per_mol === nothing
        @test_throws ArgumentError tg_activation_energy(
            :friedman, curves[1:2], targets; friedman_half_windows_min=half_windows[1:2]
        )
        @test_throws ArgumentError tg_activation_energy(:friedman, curves, targets)
        @test_throws ArgumentError tg_activation_energy(
            :starink, curves, targets; ramp_rates_K_per_min=[2.0, -5.0, 10.0]
        )
        @test_throws ArgumentError tg_activation_energy(
            :advanced_vyazovkin, curves, [0.0]; advanced=settings(0.02)
        )
    end
end
