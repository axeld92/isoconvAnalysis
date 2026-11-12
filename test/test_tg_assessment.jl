using TOML

include(joinpath(@__DIR__, "..", "scripts", "assess_tg.jl"))

@testset "TG assessment prediction baseline" begin
    prediction = TGAssessment.MassPrediction
    settings = TOML.parsefile(joinpath(@__DIR__, "..", prediction.PP_PROTOCOL))
    control = prediction.check_analytic_prediction(settings)
    @test control["held_out_max_error_mg"] < 0.001
end

@testset "TG assessment program roster" begin
    program(id, role, kind, runs; kwargs...) = Dict{String,Any}(
        "id" => id,
        "role" => role,
        "kind" => kind,
        "runs" => runs,
        "expected_repeats" => length(runs),
        (String(k) => v for (k, v) in kwargs)...,
    )
    programs = [
        program("ramp_5", "calibration", "ramp", ["r5"]; nominal_rate_K_per_min=5),
        program("ramp_10", "calibration", "ramp", ["r10"]; nominal_rate_K_per_min=10),
        program("ramp_20", "calibration", "ramp", ["r20"]; nominal_rate_K_per_min=20),
        program(
            "mass_4",
            "mass_series",
            "ramp",
            ["m1", "m2"];
            nominal_rate_K_per_min=20,
            nominal_mass_mg=4,
        ),
        program("hold_593", "holdout", "hold", ["h1", "h2", "h3"]),
    ]
    specs = [Dict("id" => id, "path" => "$id.csv") for p in programs for id in p["runs"]]
    @test isnothing(TGAssessment.validate_programs(programs, specs))

    reused = deepcopy(programs)
    push!(reused[1]["runs"], "h3")
    reused[1]["expected_repeats"] = 2
    @test_throws ArgumentError TGAssessment.validate_programs(reused, specs)

    incomplete = deepcopy(programs)
    pop!(incomplete[end]["runs"])
    @test_throws ArgumentError TGAssessment.validate_programs(incomplete, specs)

    mixed_role = deepcopy(programs)
    mixed_role[end]["kind"] = "ramp"
    mixed_role[end]["nominal_rate_K_per_min"] = 20
    @test_throws ArgumentError TGAssessment.validate_programs(mixed_role, specs)

    wrong_role = deepcopy(programs)
    wrong_role[end]["role"] = "validation"
    @test_throws ArgumentError TGAssessment.validate_programs(wrong_role, specs)

    unknown = deepcopy(programs)
    unknown[end]["runs"][1] = "absent"
    @test_throws ArgumentError TGAssessment.validate_programs(unknown, specs)

    aliased_path = deepcopy(specs)
    aliased_path[end]["path"] = aliased_path[1]["path"]
    @test_throws ArgumentError TGAssessment.validate_programs(programs, aliased_path)

    aliased_hash = deepcopy(specs)
    aliased_hash[end]["sha256"] = aliased_hash[1]["sha256"] = repeat("a", 64)
    @test_throws ArgumentError TGAssessment.validate_programs(programs, aliased_hash)
end

@testset "TG assessment from an inline run catalog" begin
    mktempdir() do directory
        audit = TGAssessment.ICTACTGAudit
        controls = audit.PublicTGEnergy.TGEnergyControls
        options = audit.load_toml(audit.PublicTGEnergy.ENERGY_PROTOCOL)
        options["controls"] = [first(options["controls"])]
        fixtures = controls.generate_tg_energy_controls(directory, options)
        config = Dict{String,Any}(
            "schema_version" => 1,
            "id" => "inline_control",
            "label" => "Known energy from inline run catalog",
            "columns" => Dict(
                "time" => "time_min",
                "time_unit" => "min",
                "temperature" => "temperature_K",
                "temperature_unit" => "K",
                "mass" => "mass_mg",
                "mass_unit" => "mg",
            ),
            "conversion" =>
                Dict("initial" => Dict("mass" => 12.0), "terminal" => Dict("mass" => 2.0)),
            "runs" => [Dict("id" => f["id"], "path" => f["filename"]) for f in fixtures],
            "programs" => [
                Dict(
                    "id" => "ramp_$i",
                    "role" => "calibration",
                    "kind" => "ramp",
                    "runs" => [f["id"]],
                    "nominal_rate_K_per_min" => f["rate_K_per_min"],
                ) for (i, f) in enumerate(fixtures)
            ],
        )
        config_path = joinpath(directory, "study.toml")
        open(io -> TOML.print(io, config), config_path, "w")
        output = joinpath(directory, "assessment")
        TGAssessment.assess_tg(config_path, output)
        report = TOML.parsefile(joinpath(output, "tg_assessment.toml"))
        @test occursin(
            "Known energy from inline run catalog",
            read(joinpath(output, "tg_assessment.md"), String),
        )
        @test !haskey(report, "prediction")
        @test length(report["programs"]) == 4
        @test length(report["repeatability"]) == 4
        for profile in report["assessment"]["baseline"]
            midpoint = only(p for p in profile["points"] if p["alpha"] == 0.50)
            @test midpoint["status"] == "estimated"
            @test isapprox(midpoint["energy_kJ_per_mol"], 160.0; atol=0.15, rtol=0)
        end

        collision = joinpath(directory, "tg_assessment.toml")
        original = read(config_path)
        write(collision, original)
        @test_throws ArgumentError TGAssessment.assess_tg(collision, directory)
        @test read(collision) == original
    end
end

# Compare existing CLI outputs when the optional public datasets have been acquired.
tg_assessment_results = get(ENV, "ISOCONV_TG_ASSESSMENT_RESULTS", "")
if !isempty(tg_assessment_results)
    @testset "Public TG assessment CLI outputs" begin
        project_root = normpath(joinpath(@__DIR__, ".."))
        pine_archive = TOML.parsefile(
            joinpath(project_root, "docs/audits/public_pine_tg.toml")
        )
        old_assessments = TOML.parsefile(
            joinpath(project_root, "docs/audits/ictac_tg_assessment.toml")
        )
        pp_archive = only(
            c for c in old_assessments["cohorts"] if c["id"] == "polypropylene"
        )
        function read_assessment(dataset)
            directory = joinpath(tg_assessment_results, dataset)
            @test occursin("ICTAC", read(joinpath(directory, "tg_assessment.md"), String))
            return TOML.parsefile(joinpath(directory, "tg_assessment.toml"))
        end
        pine = read_assessment("pine")
        pp = read_assessment("polypropylene")

        @testset "Apparent activation energy" begin
            for (report, reference) in
                ((pine, pine_archive["assessment"]), (pp, pp_archive))
                profiles = report["assessment"]["baseline"]
                @test length(profiles) == length(reference["baseline"]) == 2
                for profile in profiles
                    expected = only(
                        p for p in reference["baseline"] if p["method"] == profile["method"]
                    )
                    @test Set(profile["run_ids"]) == Set(expected["run_ids"])
                    @test length(profile["points"]) == length(expected["points"]) == 19
                    for (actual, prior) in zip(profile["points"], expected["points"])
                        @test actual["alpha"] == prior["alpha"]
                        @test actual["status"] == prior["status"]
                        if prior["status"] == "estimated"
                            @test isapprox(
                                actual["energy_kJ_per_mol"],
                                prior["energy_kJ_per_mol"];
                                atol=1e-5,
                                rtol=0,
                            )
                        end
                    end
                end
            end
            pp_checks = Dict(c["id"] => c["status"] for c in pp["assessment"]["checks"])
            @test pp_checks["repeatability"] == "missing"
            @test pp_checks["calibration"] == "missing"
            @test pp_checks["conversion"] == "partial"
        end

        @testset "Pine 593 K prediction and disjoint run roles" begin
            prediction = pine["prediction"]
            expected = pine_archive["prediction"]
            hold_ids = Set("TUBS_Wood_TGA_N2_593Kiso_R$i" for i in 1:3)
            train_ids = Set(prediction["training_ids"])
            @test length(prediction["training_ids"]) == length(train_ids) == 18
            @test train_ids == Set(expected["training_ids"])
            @test length(prediction["held_out_ids"]) == 3
            @test Set(prediction["held_out_ids"]) == hold_ids
            @test isdisjoint(train_ids, hold_ids)
            @test Set(r["run_id"] for r in pine["assessment"]["runs"]) == train_ids
            groups = pine["programs"]
            held_program = only(p for p in groups if p["role"] == "holdout")
            @test Set(held_program["run_ids"]) == hold_ids
            mass_series_ids = Set(
                id for p in groups if p["role"] == "mass_series" for id in p["run_ids"]
            )
            @test length(mass_series_ids) == 12
            @test isdisjoint(train_ids, mass_series_ids)

            scores = prediction["hold_scores"]
            @test length(scores) == 3
            @test Set(r["held_out_run"] for r in scores) == hold_ids
            parameters = (
                "energy_at_half_conversion_kJ_per_mol",
                "log_reference_rate_min_inv",
                "reaction_order",
                "retained_fraction",
            )
            @test all(r -> all(k -> r[k] == first(scores)[k], parameters), scores)
            # Refitting can pick a slightly different optimum on a newer runtime.
            # The cross-version mass budget is one increment of the source balance.
            refit_mass_atol = VERSION >= v"1.12" ? 1e-4 : 1e-7
            config_dir = joinpath(project_root, "config")
            config = TOML.parsefile(joinpath(config_dir, "tg_assessment_pine.toml"))
            specs, _ = TGAssessment.read_tg_roster(config, config_dir)
            runs = TGAssessment.import_tg_runs(config, config_dir, specs)
            policy = config["prediction"]
            mass_prediction = TGAssessment.MassPrediction
            for actual in scores
                prior = only(
                    r for r in expected["hold_scores"] if
                    r["held_out_run"] == actual["held_out_run"]
                )
                @test actual["program_id"] == held_program["id"]
                @test Set(actual["training_run_ids"]) == train_ids
                @test actual["observed_final_mass_mg"] == prior["observed_final_mass_mg"]
                @test isapprox(
                    actual["predicted_final_mass_mg"],
                    prior["predicted_final_mass_mg"];
                    atol=refit_mass_atol,
                    rtol=0,
                )
                # Replay the archived fit to test prediction independently of optimizer drift.
                held = TGAssessment.prediction_input(
                    runs[actual["held_out_run"]], policy, policy["holdout_window"]
                )
                fit = (;
                    parameters=[prior[key] for key in parameters[1:3]],
                    residue=prior["retained_fraction"],
                    shape=nothing,
                )
                replayed = mass_prediction.predict_pp(
                    fit, held.program, expected["settings"]; refined=true
                )
                @test isapprox(
                    last(replayed), prior["predicted_final_mass_mg"]; atol=1e-7, rtol=0
                )
                time = held.program.history.time_min
                last_hour =
                    searchsortedfirst(time, last(time) - policy["final_window_min"]):length(
                        time
                    )
                windows = (
                    ("evaluation", "hold", held.reaction_rows),
                    ("full", "full", eachindex(replayed)),
                    ("final_window", "final_hour", last_hour),
                )
                for (window, old_window, rows) in windows
                    @test actual[window]["rows"] == prior[old_window]["rows"]
                    replay_score = mass_prediction.score_pp(replayed, held, rows)
                    @test replay_score["rows"] == prior[old_window]["rows"]
                    for (metric, scale) in (
                        ("rmse_percent_first_mass", 100 / prior["first_mass_mg"]),
                        ("rmse_mg", 1),
                        ("mae_mg", 1),
                        ("bias_mg", 1),
                    )
                        @test isapprox(
                            actual[window][metric],
                            prior[old_window][metric];
                            atol=refit_mass_atol * scale,
                            rtol=0,
                        )
                        @test isapprox(
                            replay_score[metric],
                            prior[old_window][metric];
                            atol=1e-7 * scale,
                            rtol=0,
                        )
                    end
                end
            end
        end
    end
end
