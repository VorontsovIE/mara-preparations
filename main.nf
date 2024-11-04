#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // Stage 0: Download Data
    download_data()

    // Stage 1: Preprocess Data
    preprocessed_bed_file = preprocess_data()

    // Stage 2: Generate Flanking Regions
    flanks_file = preprocessed_bed_file
        | flanking_regions

    // Stage 3: Extract FASTA Sequences
    fasta_file = flanks_file
        | extract_fasta

    // Stage 4: Motif Analysis
    motif_result = fasta_file
        | motif_analysis
}

// Process Definitions

process download_data {
    container 'mara-preprocessing-app:latest'
    tag 'download_data'

    output:
    file 'stages/stage_00/*'

    script:
    """
    python app/MaraPreprocessingApp.py --config ${params.config_file} --stage download_data
    """
}

process preprocess_data {
    container 'mara-preprocessing-app:latest'
    tag 'preprocess_data'

    output:
    file 'stages/stage_01/hg38_promoters_v1.bed' into preprocessed_bed_file

    script:
    """
    python app/MaraPreprocessingApp.py --config ${params.config_file} --stage preprocess
    """
}

process flanking_regions {
    container 'mara-preprocessing-app:latest'
    tag 'flanking_regions'
    cpus 2

    input:
    file bed_file from preprocessed_bed_file

    output:
    file 'stages/stage_02/flanks_flanks.bed' into flanks_file

    script:
    """
    python app/MaraPreprocessingApp.py --config ${params.config_file} \
        --stage flanking_regions \
        --flank_5 ${params.flank_5} \
        --flank_3 ${params.flank_3}
    """
}

process extract_fasta {
    container 'mara-preprocessing-app:latest'
    tag 'extract_fasta'
    cpus 2

    input:
    file flanks_bed from flanks_file

    output:
    file 'stages/stage_03/sequences.fa' into fasta_file

    script:
    """
    python app/MaraPreprocessingApp.py --config ${params.config_file} --stage extract_fasta
    """
}"""
}

process motif_analysis {
    container 'mara-preprocessing-app:latest'
    tag 'motif_analysis'
    cpus 4

    input:
    file fasta_file from fasta_file

    output:
    file 'stages/stage_04/motif_occupancies.tsv' into motif_result

    script:
    """
    python app/MaraPreprocessingApp.py --config ${params.config_file} --stage motif_analysis
    """
}