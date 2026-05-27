process PREPARE_DREP_INPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    // Per-sample renamer + genomeInfo slice for cohort dereplication.
    //
    // Binette emits bins as `binette_bin{N}.fa` with no sample prefix, so
    // dropping every sample's bins into one directory for dRep would collide
    // (every sample has a `binette_bin1.fa`). Rename to `<sample>.<orig>.fa`
    // so filenames are globally unique AND traceable back to the source sample.
    //
    // Also emits the per-sample slice of the dRep --genomeInfo CSV
    // (genome,completeness,contamination) with the renamed filenames so the
    // `genome` column matches what dRep sees in input_fastas/.

    input:
    tuple val(meta), path(bins_dir, stageAs: 'bins_in/'), path(qc_tsv)

    output:
    tuple val(meta), path("renamed/*.fa"),          emit: bins
    tuple val(meta), path("${meta.id}_ginfo.csv"),  emit: genome_info

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p renamed

    # Rename each bin to <sample>.<orig>.fa. Copy (not symlink) so the staged
    # outputs are portable when publishDir uses mode 'copy' to a separate FS.
    shopt -s nullglob
    for f in bins_in/*.fa; do
        bn=\$(basename "\$f")
        cp "\$f" "renamed/${meta.id}.\${bn}"
    done

    # genomeInfo CSV slice. Binette quality TSV columns:
    #   1=name 2=origin 3=is_original 4=original_name 5=completeness 6=contamination ...
    # dRep wants: genome,completeness,contamination
    awk -F'\\t' -v OFS=',' -v sample="${meta.id}" '
        NR==1 { print "genome,completeness,contamination"; next }
        { print sample "." \$1 ".fa", \$5, \$6 }
    ' ${qc_tsv} > ${meta.id}_ginfo.csv
    """
}
