rule normalize_data:
    input:
        "Input/{file}/data.csv"
    output:
        "Output/{file}/data_normalized.csv",
        "Output/{file}/1D_penalty.png"
    script:
        "Pipeline/normalize_data.py"

rule generate_matrix:
    input:
        "Output/{file}/data_normalized.csv"
    output:
        "Output/{file}/2D_penalty.png"
    script:
        "Pipeline/generate_matrix.py"

rule compare_different_thresholds:
    input:
        "Output/{file}/data_normalized.csv"
    output:
        "Output/{file}/tolerance.png",
        "Output/{file}/mismatches.png"
    script:
        "Pipeline/compare_different_thresholds.py"