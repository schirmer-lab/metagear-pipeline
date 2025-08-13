process MSPMINER_MSPMINER{
    tag "$meta.id"
    // to review: (1) label? (2) mspminer version hardcoded
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/25.06.05/mspminer_1.1.3.sif':
        'docker.io/schirmerlab/mspminer:1.1.3' }"

    input:
        tuple val(meta), path(count_fp)

    output:
        tuple val(meta), path("mspminer"), emit: mspminer_result
        tuple val(meta), path("mspminer/all_msps.tsv"), emit: mspminer_main_table
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: ''
        def mspminer_result_dir = "mspminer"
        """

        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1

        echo RUNNING MSPMINER $meta.id
        mkdir -p ${mspminer_result_dir}
        echo "[input]
        count_matrix_file=${count_fp}
        count_matrix_has_header=true

        [genes_comparison]
        detection_limit=6
        min_prevalence=3
        min_prevalence_for_outliers=5
        max_outliers_fraction=0.3

        [genes_bins_creation]
        enabled=true
        normalize_counts=true

        [seeds_creation]
        num_best_representative_genes=30
        min_concordance=0.8
        min_seed_size=50

        [seeds_merging]
        min_concordance=0.80
        min_seed_size=150

        [msp_creation]
        min_concordance=0.8
        min_robust_concordance=0.9

        [output]
        output_dir=${mspminer_result_dir}
        print_genes_bins=true
        print_seeds=true" > ./init.txt

        mspminer ./init.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mspminer: 1.1.3
        END_VERSIONS
        """
}
