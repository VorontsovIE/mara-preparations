nextflow.enable.dsl=2

motifs_ch = Channel.fromPath("${params.motif_dir}/pwm/*.pwm").map { file -> file.baseName }
source_data_ch = Channel.value(file("source_data/"))
scripts_ch = Channel.value(file("src/"))

process download_data {
    container 'mara-preprocessing-app'
    tag 'download_data'

    input:
    path source_data

    output:
    path 'source_data/*', emit: downloaded_files // TODO Make changeable

    script:
    """
    python $params.script_file --config $params.config_file --stage download_data
    """
}

process preprocess_data {
    container 'mara-preprocessing-app'
    tag 'preprocess_data'

    input:
    path downloaded_files

    output:
    path 'stages/', emit: preprocessed_dir // TODO Simplify directory structure

    script:
    """
    python $params.script_file --config $params.config_file --tss_clusters_file $params.tss_clusters_file --stage preprocess
    """
}

process flanking_regions {
    container 'mara-preprocessing-app'
    tag 'flanking_regions'

    input:
    path preprocessed_dir

    output:
    path 'stages/', emit: flanks_bed

    script:
    """
    python $params.script_file --config $params.config_file \
        --stage flanking_regions \
        --flank_5 ${params.flank_5} \
        --flank_3 ${params.flank_3}
    """
}

process extract_fasta {
    container 'mara-preprocessing-app'
    tag 'extract_fasta'

    input:
    path genome_file // from '/work/source_data/genome/hg38.fa'
    path flanks_bed

    output:
    path 'stages/', emit: fasta_files

    script:
    """
    ls /work
    python $params.script_file --config $params.config_file --genome_file $params.genome_file --stage extract_fasta
    """
}

process motif_analysis {
    container 'mara-preprocessing-app'
    tag { "${motif}" }
    cpus 4 // TODO Make changeable

    input:
    val motif
    path scripts
    path source_data
    path fasta_files

    output:
    path "stages/stage_04/${motif}_occupancies.fasta", emit: motif_fasta
    path "stages/stage_04/${motif}_occupancies.tsv", emit: motif_tsv

    script:
    """
    python $params.script_file --config $params.config_file \
        --stage motif_analysis \
        --motif ${motif} \
    """
}

process combine_results {
    tag 'combine_results'

    input:
    path motif_fasta_list
    path motif_tsv_list

    output:
    path 'combined_motif_results.fasta'
    path 'combined_motif_results.tsv'

    script:
    """
    cat ${motif_fasta_list.flatten().join(' ')} > combined_motif_results.fasta
    cat ${motif_tsv_list.flatten().join(' ')} > combined_motif_results.tsv
    """
}

workflow {
    download_data(source_data_ch)
    preprocess_data(download_data.out.downloaded_files)
    flanking_regions(preprocess_data.out.preprocessed_dir)
    extract_fasta(source_data_ch, flanking_regions.out.flanks_bed)
    motif_analysis(motifs_ch, scripts_ch, source_data_ch, extract_fasta.out.fasta_files)
    motif_fasta_list = motif_analysis.out.motif_fasta.collect()
    motif_tsv_list = motif_analysis.out.motif_tsv.collect()
    combine_results(motif_fasta_list, motif_tsv_list)
}
