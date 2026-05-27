process GATHER_GENOMES_DIR {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    // GTDBTK_CLASSIFYWF takes a *directory* of genomes (--genome_dir). dRep
    // emits its cluster winners via `path("dereplicated_genomes/*")`, which
    // Nextflow surfaces as a list of files, not a directory. This process
    // gathers that list back into a single directory so we can hand it to
    // GTDB-Tk as a path.

    input:
    tuple val(meta), path(genomes, stageAs: 'in_genomes/*')

    output:
    tuple val(meta), path('genome_dir/'), emit: dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p genome_dir
    # `in_genomes/` is the staged copies of dRep's representatives.
    shopt -s nullglob
    for f in in_genomes/*; do
        # cp -L follows any symlinks Nextflow set up at staging so the gathered
        # dir is self-contained.
        cp -L "\$f" genome_dir/
    done
    """
}
