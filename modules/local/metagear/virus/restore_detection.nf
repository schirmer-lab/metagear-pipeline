process RESTORE_VIRAL_DETECTION {
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/25.08.11/python_3.10.sif':
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    path per_sample_dir // a previous run's virus/per_sample

    output:
    path "*_Merged_Genomad_CheckV_Summary.tsv"          , emit: virus_tables             , optional: true
    path "*_Merged_Genomad_CheckV_Summary.filtered.tsv" , emit: virus_filtered_tables    , optional: true
    path "*_Plasmid_Summary.filtered.tsv"               , emit: plasmid_filtered_tables  , optional: true
    path "*_viral_ids_to_keep.txt"                      , emit: viral_ids                , optional: true
    path "*_plasmid_ids_to_keep.txt"                    , emit: plasmid_ids              , optional: true
    path "*.virus.filtered.gz"                          , emit: viral_sequences          , optional: true
    path "*.plasmid.filtered.gz"                        , emit: plasmid_sequences        , optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Undo the rename that publishing applied, so every consumer downstream sees the
    # names VIRAL_DETECTION itself emitted. Published names collide across samples
    # (the sample sits in the directory, not the file name) and virus.filtered.fna.gz
    # also matches the output glob of the process that concatenates it.
    restore() {
        if [ -f "\$1" ]; then
            cp -L "\$1" "\$2"
        fi
    }

    for sample_dir in ${per_sample_dir}/*/ ; do
        sample=\$(basename "\$sample_dir")

        restore "\$sample_dir/virus.summary.tsv"            "\${sample}_Merged_Genomad_CheckV_Summary.tsv"
        restore "\$sample_dir/virus.summary.filtered.tsv"   "\${sample}_Merged_Genomad_CheckV_Summary.filtered.tsv"
        restore "\$sample_dir/plasmid.summary.filtered.tsv" "\${sample}_Plasmid_Summary.filtered.tsv"
        restore "\$sample_dir/virus.ids.txt"                "\${sample}_viral_ids_to_keep.txt"
        restore "\$sample_dir/plasmid.ids.txt"              "\${sample}_plasmid_ids_to_keep.txt"
        restore "\$sample_dir/virus.filtered.fna.gz"        "\${sample}.virus.filtered.gz"
        restore "\$sample_dir/plasmid.filtered.fna.gz"      "\${sample}.plasmid.filtered.gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(cp --version 2>&1 | sed -n '1s/^cp (GNU coreutils) //p')
    END_VERSIONS
    """
}
