process STAGE_DREP_WORK {
    tag "cohort"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    // Concatenates per-sample genomeInfo CSVs (from PREPARE_DREP_INPUTS) into
    // a single cohort-level CSV and stages it inside a `drep_work_seed/`
    // directory. The DREP_DEREPLICATE module accepts this directory through
    // its second input slot (`tuple val(meta2), path(drep_work, stageAs:
    // 'drep_work/')`), so the cohort CSV lands at `drep_work/genomeInfo.csv`
    // inside dRep's work directory. ext.args on DREP_DEREPLICATE then points
    // `--genomeInfo drep_work/genomeInfo.csv`, letting dRep skip its internal
    // CheckM run.

    input:
    path csvs   // collected list of per-sample <sample>_ginfo.csv

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
    # Take the header from the first file; concat data rows from all.
    csvs_arr=( ${csvs} )
    head -n1 "\${csvs_arr[0]}" > drep_work_seed/genomeInfo.csv
    for f in ${csvs}; do
        tail -n +2 "\$f" >> drep_work_seed/genomeInfo.csv
    done
    """
}
