@testset "Observed TG conversion support" begin
    mktempdir() do directory
        function tg_curve(
            name,
            alpha;
            temperature=collect(range(300.0, 400.0; length=max(length(alpha), 2)))[1:length(
                alpha
            )],
            clock=collect(7:(6 + length(alpha))),
            material="polymer",
            unit=:mg,
            id=name,
        )
            source = joinpath(directory, "$name.csv")
            timed = clock !== nothing
            open(source, "w") do io
                println(io, timed ? "T,Time,Mass" : "T,Mass")
                println(io, timed ? "270,0,8" : "270,8")
                for i in eachindex(alpha)
                    i == 2 && println(io)
                    mass = 8 * (1 - alpha[i])
                    println(
                        io,
                        if timed
                            "$(temperature[i]),$(clock[i]),$mass"
                        else
                            "$(temperature[i]),$mass"
                        end,
                    )
                end
            end
            mass = ThermalSignalSpec(:m, ThermalColumnSpec("Mass", unit), :mass)
            mapping = ThermalColumnMapping(
                ThermalColumnSpec("T", :K),
                [mass];
                time=timed ? ThermalColumnSpec("Time", :min) : nothing,
            )
            run = read_thermal_delimited(source, mapping; id=id, material_id=material)
            return tg_conversion(
                run,
                :m;
                segment=1,
                reaction_rows=2:(length(alpha) + 1),
                initial_reference=8,
                terminal_reference=0,
            )
        end

        @testset "Crossings and exact encounters" begin
            curve = tg_curve(
                "ascending", [0, 0.75, 1]; temperature=[320, 300, 360], clock=[8, 10, 13]
            )
            bracket = tg_conversion_crossings(curve, 0.5)
            @test bracket.alpha == 0.5
            @test bracket.status === :unique_upward
            event = only(bracket.crossings)
            @test event.kind === :upward
            @test event.segment_rows == 2:3
            @test event.source_lines == (3, 5)
            @test event.temperature_range_K == (300.0, 320.0)
            @test event.temperature_K ≈ 320 - 20 * 2 / 3
            @test event.time_min ≈ 8 + 2 * 2 / 3

            exact = tg_conversion_crossings(curve, 0.75)
            @test exact.status === :unique_upward
            @test only(exact.crossings).segment_rows == 3:3
            @test only(exact.crossings).source_lines == (5, 5)
            @test only(exact.crossings).temperature_K == 300
            @test only(exact.crossings).time_min == 10
            for (target, line) in ((0.0, 3), (1.0, 6))
                endpoint = tg_conversion_crossings(curve, target)
                @test endpoint.status === :endpoint
                @test only(endpoint.crossings).source_lines == (line, line)
            end

            plateau_curve = tg_curve(
                "plateau", [0, 0.5, 0.5, 0.5, 1]; temperature=[300, 340, 390, 320, 450]
            )
            plateau = tg_conversion_crossings(plateau_curve, 0.5)
            @test plateau.status === :plateau
            @test only(plateau.crossings).segment_rows == 3:5
            @test only(plateau.crossings).temperature_range_K == (320.0, 390.0)
            @test only(plateau.crossings).temperature_K === nothing
            @test only(plateau.crossings).time_min === nothing

            tangent = tg_conversion_crossings(tg_curve("tangent", [0, 0.5, 0.25]), 0.5)
            @test tangent.status === :tangency
            @test only(tangent.crossings).temperature_K == 350
            @test tg_conversion_crossings(tg_curve("single", [0.5]), 0.5).status ===
                :endpoint
            @test tg_conversion_crossings(tg_curve("all_equal", [0.5, 0.5]), 0.5).status ===
                :plateau
            down = tg_curve("downward", [0.75, 0.5, 0.25]; clock=nothing)
            @test tg_conversion_crossings(down, 0.5).status === :unique_downward
            @test tg_conversion_crossings(down, 0.625).status === :unique_downward
            @test only(tg_conversion_crossings(down, 0.625).crossings).time_min === nothing
            missing = tg_conversion_crossings(down, 0.8)
            @test missing.status === :missing && isempty(missing.crossings)

            reverse = tg_curve("reverse", [0, 0.75, 0.25, 1])
            ambiguous = tg_conversion_crossings(reverse, 0.5)
            @test ambiguous.status === :ambiguous
            @test getfield.(ambiguous.crossings, :kind) == [:upward, :downward, :upward]
            @test getfield.(ambiguous.crossings, :source_lines) == [(3, 5), (5, 6), (6, 7)]
        end

        @testset "Common coverage, basis and clock provenance" begin
            first_run = tg_curve("first", [0, 0.5, 1]; temperature=[300, 350, 400])
            second_run = tg_curve(
                "second",
                [0.25, 0.75, 1];
                temperature=[350, 450, 500],
                clock=nothing,
                unit=:g,
            )
            targets = [0.125, 0.5, 0.875]
            support = tg_conversion_support([first_run, second_run], targets)
            @test support.material_id == "polymer"
            @test support.mass_basis === :absolute_mass
            @test support.run_ids == ["first", "second"]
            @test support.missing_time_run_ids == ["second"]
            @test getfield.(support.points, :common_upward) == [false, true, true]
            unavailable, common = support.points[1:2]
            @test unavailable.upward_run_count == 1
            @test unavailable.temperature_min_K === nothing
            @test unavailable.temperature_mean_K === nothing
            @test unavailable.temperature_span_K === nothing
            @test unavailable.reciprocal_temperature_span_K_inv === nothing
            @test unavailable.reciprocal_temperature_sxx_K_inv2 === nothing
            @test common.upward_run_count == 2
            @test common.temperature_min_K == 350
            @test common.temperature_max_K == 400
            @test common.temperature_mean_K == 375
            @test common.temperature_span_K == 50
            reciprocal_span = 1 / 350 - 1 / 400
            @test common.reciprocal_temperature_span_K_inv ≈ reciprocal_span
            @test common.reciprocal_temperature_sxx_K_inv2 ≈ reciprocal_span^2 / 2
            @test [only(d.crossings).temperature_K for d in common.crossings] == [350, 400]

            duplicate = tg_curve("duplicate", [0, 0.5, 1]; temperature=[300, 350, 400])
            same_T = tg_curve(
                "same_T", [0, 0.5, 1]; temperature=[300, 350, 400], clock=[1, 2, 3]
            )
            @test only(tg_conversion_support([first_run, same_T], [0.5]).points).temperature_span_K ==
                0
            @test_throws ArgumentError tg_conversion_support(
                [first_run, duplicate], targets
            )
            same_id = tg_curve(
                "same_id", [0, 0.5, 1]; temperature=[300, 400, 420], id="first"
            )
            @test_throws ArgumentError tg_conversion_support([first_run, same_id], targets)

            unknown = tg_curve("unknown", [0, 1]; material=nothing)
            other = tg_curve("other", [0, 1]; material="other")
            relative = tg_curve("relative", [0, 1]; unit=:percent)
            @test_throws ArgumentError tg_conversion_support([first_run, unknown], targets)
            @test tg_conversion_support(
                [first_run, unknown], targets; material_id="polymer"
            ).material_id == "polymer"
            @test_throws ArgumentError tg_conversion_support(
                [first_run, other], targets; material_id="polymer"
            )
            @test_throws ArgumentError tg_conversion_support([first_run, relative], targets)
            @test_throws ArgumentError tg_conversion_support(TGConversionResult[], targets)
            @test_throws ArgumentError tg_conversion_support([first_run], Float64[])
            @test_throws ArgumentError tg_conversion_support([first_run], [0.5, 0.5])
            @test_throws ArgumentError tg_conversion_support([first_run], [0.75, 0.5])
            @test_throws ArgumentError tg_conversion_crossings(first_run, NaN)
            @test_throws ArgumentError tg_conversion_crossings(first_run, 1.1)
            first_run.alpha[2] = NaN
            @test_throws ArgumentError tg_conversion_crossings(first_run, 0.5)
        end
    end
end
