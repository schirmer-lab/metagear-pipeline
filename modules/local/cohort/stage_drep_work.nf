process STAGE_DREP_WORK {
    tag "cohort"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    // Builds the cohort-level genomeInfo CSV that dRep consumes via
    // --genomeInfo. Reads per-sample Binette QC TSVs directly (eliminating
    // the now-redundant PREPARE_DREP_INPUTS per-sample CSV slice step that
    // previously sat between BINETTE and here) and transforms in one pass:
    //
    //   Binette QC TSV cols: 1=name 2=origin 3=is_original 4=original_name
    //                       5=completeness 6=contamination 7=score ...
    //   dRep genomeInfo CSV: genome,completeness,contamination
    //
    // Binette is invoked with `--prefix <sample>` upstream, so col 1 already
    // reads `<sample>_binN` (matches the FASTA filename minus `.fa`). The
    // awk just needs to append `.fa` and reorder.
    //
    // The CSV lands at `drep_work_seed/genomeInfo.csv` and gets staged inside
    // dRep's work directory through DREP_DEREPLICATE's second input slot
    // (`stageAs: 'drep_work/'`). ext.args on DREP_DEREPLICATE then points
    // `--genomeInfo drep_work/genomeInfo.csv`, letting dRep skip CheckM.

    input:
    // Per-sample QC TSVs are published as `<sample>.quality_report.tsv` so
    // filenames are globally unique and self-descriptive when shared. They
    // stage cleanly here as a flat list.
    path qc_tsvs

    output:
    // Plain path — Nextflow's parser refuses map literals inside val() at
    // module level ("No such variable: id"). The consumer subworkflow wraps
    // this in a meta tuple via .map { dir -> [[id: 'drep_work'], dir] } before
    // passing it to DREP_DEREPLICATE.
    path('drep_work_seed/'), emit: drep_work

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p drep_work_seed
    awk -F'\\t' -v OFS=',' '
        BEGIN { print "genome,completeness,contamination" }
        FNR==1 { next }
        { print \$1 ".fa", \$5, \$6 }
    ' ${qc_tsvs} > drep_work_seed/genomeInfo.csv
    """
}
