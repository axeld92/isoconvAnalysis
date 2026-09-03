using MAT

const REPRESENTATIVE_DYNAMIC_MATRIX = [
    25.0 0.0 0.10 100.0 80.0 20.0 0.92
    100.0 1.0 0.20 98.0 80.0 20.0 0.90
    200.0 2.0 0.30 90.0 80.0 20.0 0.85
    NaN NaN NaN NaN NaN NaN NaN
]

const REPRESENTATIVE_RAMP_HOLD_MATRIX = [
    25.0 0.0 0.10 100.0 20.0 20.0 0.92 1.0
    100.0 1.0 0.20 98.0 20.0 20.0 0.90 3.0
    120.0 2.0 0.30 90.0 20.0 20.0 0.85 3.0
    120.0 3.0 0.25 88.0 20.0 20.0 0.85 4.0
    NaN NaN NaN NaN NaN NaN NaN NaN
]

function write_representative_mat(path::AbstractString)
    MAT.matwrite(
        path,
        Dict(
            "dynamic_fixture" => REPRESENTATIVE_DYNAMIC_MATRIX,
            "ramp_hold_fixture" => REPRESENTATIVE_RAMP_HOLD_MATRIX,
            "not_a_matrix" => "fixture metadata",
        ),
    )
    return path
end
