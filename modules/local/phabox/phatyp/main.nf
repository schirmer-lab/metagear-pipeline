process PHABOX2_PHATYP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/phabox:2.1.13--pyhdfd78af_1':
        'quay.io/biocontainers/phabox:2.1.13--pyhdfd78af_1' }"

    // Lifestyle prediction (temperate vs virulent) for viral contigs.
    //
    // Scoped deliberately to `--task phatyp`. PhaBOX2 also ships `cherry` (host
    // prediction), `phagcn` (taxonomy) and `phavip` (protein annotation), which
    // would duplicate IPHOP_PREDICT, GENOMAD and PHAROKKA/PHOLD respectively.
    // Running `end_to_end` here would produce three competing annotations for
    // quantities this pipeline already derives elsewhere.
    //
    // PhaTYP is a BERT model over protein-cluster "sentences", so unlike BACPHLIP
    // it is built for fragmentary contigs rather than complete genomes. That is
    // what makes it usable on an assembly-derived catalog whose median
    // representative is well under 2 kb; BACPHLIP's own documentation states it
    // must not be run on fragmented genomes, which is why it is not wrapped here.

    input:
    tuple val(meta), path(fasta)
    path phatyp_db

    output:
    tuple val(meta), path("*.phatyp_prediction.tsv"), emit: lifestyle
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # phabox2 reads plain FASTA; the catalog chunks arrive gzipped from
    # SEQKIT_SPLIT2.
    INPUT=$fasta
    if [[ $fasta == *.gz ]]; then
        gunzip -c $fasta >| \$PWD/${prefix}_plain.fasta
        INPUT=\$PWD/${prefix}_plain.fasta
    fi

    phabox2 \\
        --task phatyp \\
        --contigs \$INPUT \\
        --outpth phatyp_results \\
        --dbdir $phatyp_db \\
        --threads $task.cpus \\
        $args

    # phabox2 writes final_prediction/phatyp_prediction.tsv (columns Accession,
    # Length, TYPE, PhaTYPScore) plus a phavip_prediction.tsv produced as a
    # by-product of its protein step. Only the lifestyle table is emitted; the
    # protein annotation duplicates PHAROKKA/PHOLD output.
    #
    # TYPE carries four values, not two, and the distinction matters downstream:
    # `temperate` and `virulent` are calls, `filtered` marks a contig below the
    # length threshold that was never classified, and `-` marks one that was
    # attempted and returned no call. Only the last two are non-calls and only
    # `filtered` should leave the denominator of a temperate fraction.
    if [[ -s phatyp_results/final_prediction/phatyp_prediction.tsv ]]; then
        mv phatyp_results/final_prediction/phatyp_prediction.tsv ./${prefix}.phatyp_prediction.tsv
    else
        # A chunk whose contigs all fall below the length threshold produces no
        # table. Emit the header alone, so the merge sees a well-formed empty
        # result and a genuinely failed run stays distinguishable from a
        # legitimately empty one.
        printf 'Accession\\tLength\\tTYPE\\tPhaTYPScore\\n' > ./${prefix}.phatyp_prediction.tsv
    fi

    # phabox2 has no --version flag. The version appears in the banner every
    # invocation prints, so it is read from --help; note the banner is absent from
    # the task-specific help (`--task phatyp -h`), which prints only the flag list.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phabox: \$(phabox2 --help 2>&1 | sed -e 's/\\x1b\\[[0-9;]*m//g' | grep -om1 'PhaBOX v[0-9.]*' | sed 's/^PhaBOX v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Accession\\tLength\\tTYPE\\tPhaTYPScore\\n'  > ${prefix}.phatyp_prediction.tsv
    printf 'stub_contig_1\\t45000\\ttemperate\\t0.98\\n' >> ${prefix}.phatyp_prediction.tsv
    printf 'stub_contig_2\\t12000\\tvirulent\\t0.91\\n'  >> ${prefix}.phatyp_prediction.tsv
    printf 'stub_contig_3\\t800\\tfiltered\\t0\\n'       >> ${prefix}.phatyp_prediction.tsv
    printf 'stub_contig_4\\t1200\\t-\\t0\\n'             >> ${prefix}.phatyp_prediction.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phabox: 2.1.13
    END_VERSIONS
    """
}
