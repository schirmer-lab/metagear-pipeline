process EXTRACT_GENES {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/python310:25.09.10' :
        'docker.io/schirmerlab/python310:25.09.10' }"

    input:
    tuple val(meta), path(contig_ids), path(all_genes)

    output:
    tuple val(meta), path("*.fasta"), emit: extracted_genes, optional: true
    tuple val(meta), path("*.ids.txt"), emit: extracted_gene_ids, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Split contig_ids into full contigs and partial (with boundaries)
    grep -v '|provirus' ${contig_ids} > full_contigs.txt || touch full_contigs.txt
    grep '|provirus' ${contig_ids} > partial_contigs.txt || touch partial_contigs.txt
    
    # Process full contigs
    if [ -s full_contigs.txt ]; then
        # Create grep patterns from full contig IDs (add ^> prefix for fasta headers)
        sed 's/^/^>/g' full_contigs.txt > full_contig_patterns.txt
        
        # Extract matching gene headers and remove > prefix for gene IDs
        zgrep -f full_contig_patterns.txt ${all_genes} | sed 's/^>//' > ${prefix}.full.txt
        
        # Use seqtk subseq to extract sequences - fast and efficient
        seqtk subseq ${all_genes} ${prefix}.full.txt > ${prefix}.full.fasta.tmp
    else
        # Create empty files if no full contigs
        touch ${prefix}.full.txt ${prefix}.full.fasta.tmp
    fi
    
    # Process partial contigs (with boundaries)
    if [ -s partial_contigs.txt ]; then
        # Step 1: Extract just the contig names (without boundaries) for fast filtering
        cut -d'|' -f1 partial_contigs.txt | sed 's/^/^>/g' > partial_contig_patterns.txt
        
        # Step 2: Fast filter all_genes to get just headers from relevant contigs
        zgrep -f partial_contig_patterns.txt ${all_genes} > filtered_headers.txt
        
        # Step 3: Apply boundary logic on the much smaller filtered dataset
        cat filtered_headers.txt | awk -F'::' -v contig_file="partial_contigs.txt" '
        BEGIN {
            # Read contig boundaries
            while ((getline line < contig_file) > 0) {
                split(line, parts, "|")
                contig = parts[1]
                split(parts[2], boundary_parts, "_")
                prov_start = boundary_parts[2]
                prov_end = boundary_parts[3]
                boundaries[contig] = prov_start ":" prov_end
            }
            close(contig_file)
        }
        {
            # Parse gene header: >P13752_344_S66_L001_k141_40705::28::31043::31828::-
            header = \$0
            gsub(/^>/, "", header)  # Remove >
            
            # Extract components
            contig = \$1
            gene_start = \$3
            gene_end = \$4
            
            # Check if this contig has boundaries defined
            if (contig in boundaries) {
                split(boundaries[contig], bounds, ":")
                prov_start = bounds[1]
                prov_end = bounds[2]
                
                # Check if gene is completely within provirus boundaries
                if (gene_start >= prov_start && gene_end <= prov_end) {
                    print header
                }
            }
        }' > ${prefix}.partial.txt
        
        # Use seqtk subseq to extract sequences - fast and efficient
        seqtk subseq ${all_genes} ${prefix}.partial.txt > ${prefix}.partial.fasta.tmp
    else
        # Create empty files if no partial contigs
        touch ${prefix}.partial.txt ${prefix}.partial.fasta.tmp
    fi
    
    # Combine all gene IDs and create combined fasta
    cat ${prefix}.full.txt ${prefix}.partial.txt > ${prefix}.ids.txt
    cat ${prefix}.full.fasta.tmp ${prefix}.partial.fasta.tmp > ${prefix}.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
