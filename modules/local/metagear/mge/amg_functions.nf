process EXTRACT_AMGS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/python310:25.09.10' :
        'docker.io/schirmerlab/python310:25.09.10' }"

    input:
    tuple val(meta), path(amg_summary), path(amg_fna), path(amg_faa)

    output:
    tuple val(meta), path("*.faa"), emit: amgs_faa, optional: true
    tuple val(meta), path("*.fna"), emit: amgs_fna, optional: true
    tuple val(meta), path("*.txt"), emit: amg_dramv_ids, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    awk -F'\\t' 'NR==1 {next} (\$4+0) <= 3 {print \$1}' ${amg_summary}| sort | uniq > amg_filtered_ids.txt

    cat ${amg_faa} | seqtk subseq - amg_filtered_ids.txt > amg_filtered.faa
    cat ${amg_fna} | seqtk subseq - amg_filtered_ids.txt > amg_filtered.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | awk '/Version:/ {print \$2; exit}')
    END_VERSIONS
    """
}


process MAP_AMG_CATALOG {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/python310:25.09.10' :
        'docker.io/schirmerlab/python310:25.09.10' }"

    input:
    tuple val(meta), path(amg_summary), path(tpm_table), path(rpkm_table), path(count_table), path(mmseqs_search_tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: amg_catalog_summary, optional: true
    tuple val(meta), path("amg_catalog_ids.txt"), emit: amg_catalog_ids, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    ## Keep the best single target from the search results
    awk -F'\\t' '
    {
    q=\$1;
    if(!(q in keep) || \$7>kb[q] || (\$7==kb[q] && \$3>kp[q])) {keep[q]=\$0; kb[q]=\$7; kp[q]=\$3}
    }
    END{
    OFS="\\t";
    for(q in keep){ split(keep[q],a,"\\t"); print a[2], a[1] }
    }' $mmseqs_search_tsv > amg_catalog_map.txt


    # Normalize newlines; set bytewise collation
    sed -i 's/\r\$//' amg_catalog_map.txt $amg_summary; export LC_ALL=C


    # Header (take AMG header minus first column)
    hdr_sub=\$(head -n1 "$amg_summary" | cut -f2-)
    printf "catalog_id\tgene\t%s\n" "\$hdr_sub" > amg_to_catalog.tsv


    # Sort on join keys (map: col2; amg: col1), skip AMG header
    sort -t \$'\\t' -k2,2 amg_catalog_map.txt > amg_catalog_map.sorted.txt
    tail -n +2 "$amg_summary" | sort -t \$'\\t' -k1,1 > "amg_summary.sorted.txt"


    # Build explicit -o field list: 2.2,2.3,...,2.N
    N=\$(head -n1 "$amg_summary" | awk -F'\t' '{print NF}')
    flds=\$(awk -v N="\$N" 'BEGIN{for(i=2;i<=N;i++){printf ",2.%d",i}}')

    # Join: keep only matches; to keep all map rows add: -a 1 -e ""
    eval "join -t \$'\\t' -1 2 -2 1 -o 1.1,1.2\${flds} \"amg_catalog_map.sorted.txt\" \"amg_summary.sorted.txt\"" >> amg_to_catalog.tsv

    cat amg_to_catalog.tsv | cut -f1 | sort | uniq > amg_catalog_ids.txt

    LC_ALL=C awk -F \$'\\t' -v OFS=\$'\\t' '
    NR==FNR { gsub(/\\r\$/,"",\$1); ids[\$1]=1; next }   # load IDs (one per line)
    FNR==1 || (\$1 in ids)                            # keep header or matching first column
    ' amg_catalog_ids.txt $tpm_table > virus.amg_tpm.tsv

    LC_ALL=C awk -F \$'\\t' -v OFS=\$'\\t' '
    NR==FNR { gsub(/\\r\$/,"",\$1); ids[\$1]=1; next }   # load IDs (one per line)
    FNR==1 || (\$1 in ids)                            # keep header or matching first column
    ' amg_catalog_ids.txt $rpkm_table > virus.amg_rpkm.tsv

    LC_ALL=C awk -F \$'\\t' -v OFS=\$'\\t' '
    NR==FNR { gsub(/\\r\$/,"",\$1); ids[\$1]=1; next }   # load IDs (one per line)
    FNR==1 || (\$1 in ids)                            # keep header or matching first column
    ' amg_catalog_ids.txt $count_table > virus.amg_count.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | awk '/Version:/ {print \$2; exit}')
    END_VERSIONS
    """
}
