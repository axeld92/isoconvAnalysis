using SHA: sha256

@testset "delimited thermal import" begin
    mass = ThermalSignalSpec(:mass, ThermalColumnSpec("mass", :mg), :mass)
    mapping = ThermalColumnMapping(
        ThermalColumnSpec("temp", :C),
        [mass];
        time=ThermalColumnSpec("time, min", :min),
        segment_column="phase",
    )
    mktempdir() do directory
        path = joinpath(directory, "thermal.run.csv")
        original =
            "\ufeff\r\ntemp,\" time, min \",mass,phase,note\r\n25,0,10,\"ramp, first\",\"a\"\"b\"\r\n\r\n" *
            "30,1,9.8,\"ramp, first\",\"\"\r\n30,0,9.7,hold,ignored\r\n27,0,9.6,\"ramp, first\",\r\n"
        write(path, original)
        fingerprint = bytes2hex(sha256(codeunits(original)))
        run = read_thermal_delimited(path, mapping; expected_sha256=uppercase(fingerprint))
        @test run.id == "thermal.run"
        @test run.material_id === nothing
        @test run.temperature == [25, 30, 30, 27]
        @test run.time == [0, 1, 0, 0]
        @test only(run.signals).values == [10, 9.8, 9.7, 9.6]
        @test getproperty.(run.segments, :label) == ["ramp, first", "hold", "ramp, first"]
        @test getproperty.(run.segments, :rows) == [1:2, 3:3, 4:4]
        @test run.provenance.source_path == abspath(path)
        @test run.provenance.source_lines == [3, 5, 6, 7]
        @test run.provenance.unmapped_columns == ["note"]
        @test run.provenance.sha256 == fingerprint
        @test_throws DataImportError read_thermal_delimited(
            path, mapping; expected_sha256=repeat("0", 64)
        )
        @test_throws DataImportError read_thermal_delimited(
            path, mapping; expected_sha256="bad"
        )

        plain_mapping = ThermalColumnMapping(ThermalColumnSpec("temp", :C), [mass])
        for contents in (
            "temp,mass,mass\n25,10,11\n",
            "temp,mass,\n25,10,11\n",
            "temp,other\n25,10\n",
            "temp,mass\n25,\n",
            "temp,mass\n25,NaN\n",
            "temp,mass\n25,Inf\n",
            "temp,mass\n25,10,\n",
            "temp,mass\n25,\"10\" extra\n",
            "temp,mass\n25,1\"0\n",
            "temp,mass\n25,\"10\n11\"\n",
            "temp,mass\n",
            "\n",
        )
            write(path, contents)
            @test_throws DataImportError read_thermal_delimited(path, plain_mapping)
        end
        write(path, "temp\tmass\n25\t10\n\t\n")
        @test_throws r"line 3: expected a finite number" read_thermal_delimited(
            path, plain_mapping; delimiter=('\t')
        )
        write(path, "temp,\"time, min\",mass,phase\n25,0,10,\" \"\n")
        @test_throws DataImportError read_thermal_delimited(path, mapping)
    end

    @testset "separate metadata separator" begin
        mktempdir() do directory
            path = joinpath(directory, "instrument.txt")
            write(path, "#DATE/TIME:19.03.2021 15:48:02\n##temp;mass\n25;10\n")
            mapping = ThermalColumnMapping(ThermalColumnSpec("temp", :C), [mass])
            run = read_thermal_delimited(
                path, mapping; delimiter=(';'), header_prefix="##", metadata_delimiter=(':')
            )
            @test run.metadata["DATE/TIME"] == "19.03.2021 15:48:02"
            @test run.provenance.source_lines == [3]
            @test only(run.signals).values == [10.0]
        end
    end

    @testset "declared column unit labels" begin
        mktempdir() do directory
            path = joinpath(directory, "units.csv")
            labels = ["[°C]", "[s]", "[mg]", "[]"]
            mapping = ThermalColumnMapping(
                ThermalColumnSpec("temp", :C),
                [mass];
                time=ThermalColumnSpec("time", :s),
                segment_column="phase",
            )
            write(path, "temp,time,mass,phase\n[°C],[s],[mg],[]\n25,60,10,ramp\n")
            run = read_thermal_delimited(path, mapping; unit_labels=labels)
            @test run.temperature == [25]
            @test time_min(run) == [1]
            @test run.provenance.source_lines == [3]
            @test_throws r"line 2: column unit labels do not match" read_thermal_delimited(
                path, mapping; unit_labels=["[K]", labels[2:end]...]
            )
            @test_throws r"line 2: expected a finite number" read_thermal_delimited(
                path, mapping
            )
            write(path, "temp,time,mass,phase\n")
            @test_throws r"missing column unit labels" read_thermal_delimited(
                path, mapping; unit_labels=labels
            )
        end
    end

    @testset "explicit Windows-1252 export" begin
        mktempdir() do directory
            path = joinpath(directory, "netzsch.csv")
            open(path, "w") do io
                write(io, "#EXO:;-1\r\n#CORR. CODE:;000\r\n#COMMENT:;\"A;B\"\r\n##Temp./")
                write(io, 0x9b)
                write(
                    io,
                    "C;TG/%;DSC/(mW/mg);Segment;Other\r\n30;99.9;-0.1;1;\r\n31;99.8;0.1;2;\r\n",
                )
            end
            dsc = ThermalSignalSpec(
                :dsc,
                ThermalColumnSpec("DSC/(mW/mg)", :mW_per_mg),
                :specific_heat_flow;
                exotherm_sign=-1,
            )
            tg = ThermalSignalSpec(:mass, ThermalColumnSpec("TG/%", :percent), :mass)
            mapping = ThermalColumnMapping(
                ThermalColumnSpec("Temp./›C", :C), [tg, dsc]; segment_column="Segment"
            )
            @test_throws DataImportError read_thermal_delimited(
                path, mapping; delimiter=(';'), header_prefix="##"
            )
            run = read_thermal_delimited(
                path,
                mapping;
                id="heat",
                material_id="sample",
                delimiter=(';'),
                header_prefix="##",
                encoding=:cp1252,
            )
            @test run.id == "heat"
            @test run.material_id == "sample"
            @test run.time === nothing
            @test run.metadata ==
                Dict("EXO" => "-1", "CORR. CODE" => "000", "COMMENT" => "A;B")
            @test getproperty.(run.segments, :label) == ["1", "2"]
            @test run.provenance.source_lines == [5, 6]
            @test run.provenance.unmapped_columns == ["Other"]
            @test run.signals[2].spec.exotherm_sign == -1
            @test run.signals[2].values == [-0.1, 0.1]
        end
    end
end
