#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <sequence_ids.txt> <genes.fasta[.gz]>" >&2
    exit 1
fi

sequence_ids=$1
genes_file=$2
prefix=${PREFIX:-extracted_genes}

if [[ ! -f "$sequence_ids" ]]; then
    echo "Sequence ID file '$sequence_ids' not found" >&2
    exit 1
fi

if [[ ! -f "$genes_file" ]]; then
    echo "Genes FASTA file '$genes_file' not found" >&2
    exit 1
fi

full_ids_file="${prefix}.full.txt"
partial_ids_file="${prefix}.partial.txt"
full_fasta_tmp="${prefix}.full.fasta.tmp"
partial_fasta_tmp="${prefix}.partial.fasta.tmp"
all_ids_file="${prefix}.ids.txt"
combined_fasta="${prefix}.fasta"
filtered_headers="filtered_headers.txt"
full_contigs="full_contigs.txt"
partial_contigs="partial_contigs.txt"
full_patterns="full_contig_patterns.txt"
partial_patterns="partial_contig_patterns.txt"

> "$full_contigs"
> "$partial_contigs"

if ! grep -v '|provirus' "$sequence_ids" > "$full_contigs"; then
    :
fi
if ! grep '|provirus' "$sequence_ids" > "$partial_contigs"; then
    :
fi

if [[ "$genes_file" == *.gz ]]; then
    grep_cmd=(zgrep -E)
else
    grep_cmd=(grep -E)
fi

if [[ -s "$full_contigs" ]]; then
    awk '{ printf("^>%s(::|$)\n", $0) }' "$full_contigs" > "$full_patterns"
    if ! "${grep_cmd[@]}" -f "$full_patterns" "$genes_file" | sed 's/^>//' > "$full_ids_file"; then
        : > "$full_ids_file"
    fi
    if [[ -s "$full_ids_file" ]]; then
        seqtk subseq "$genes_file" "$full_ids_file" > "$full_fasta_tmp"
    else
        > "$full_fasta_tmp"
    fi
else
    > "$full_patterns"
    > "$full_ids_file"
    > "$full_fasta_tmp"
fi

if [[ -s "$partial_contigs" ]]; then
    cut -d'|' -f1 "$partial_contigs" | awk '{ printf("^>%s(::|$)\n", $0) }' > "$partial_patterns"
    if ! "${grep_cmd[@]}" -f "$partial_patterns" "$genes_file" > "$filtered_headers"; then
        : > "$filtered_headers"
    fi

    awk -F'::' -v contig_file="$partial_contigs" '
        BEGIN {
            while ((getline line < contig_file) > 0) {
                split(line, parts, "|")
                if (length(parts) < 2) {
                    continue
                }
                contig = parts[1]
                split(parts[2], boundary_parts, "_")
                if (length(boundary_parts) < 3) {
                    continue
                }
                prov_start = boundary_parts[2]
                prov_end = boundary_parts[3]
                boundaries[contig] = prov_start ":" prov_end
            }
            close(contig_file)
        }
        {
            header = $0
            gsub(/^>/, "", header)
            # Strip the FASTA ">" from the contig field too. gsub above rewrites
            # a copy of the record, so $1 still carries the leading ">" and the
            # boundaries[] lookup below would never match, silently dropping
            # every provirus gene.
            contig = $1
            sub(/^>/, "", contig)
            gene_start = $3 + 0
            gene_end = $4 + 0
            if (!(contig in boundaries)) {
                next
            }
            split(boundaries[contig], bounds, ":")
            prov_start = bounds[1] + 0
            prov_end = bounds[2] + 0
            if (gene_start >= prov_start && gene_end <= prov_end) {
                print header
            }
        }
    ' "$filtered_headers" > "$partial_ids_file"

    if [[ -s "$partial_ids_file" ]]; then
        seqtk subseq "$genes_file" "$partial_ids_file" > "$partial_fasta_tmp"
    else
        > "$partial_fasta_tmp"
    fi
else
    > "$partial_patterns"
    > "$filtered_headers"
    > "$partial_ids_file"
    > "$partial_fasta_tmp"
fi

cat "$full_ids_file" "$partial_ids_file" > "$all_ids_file"
cat "$full_fasta_tmp" "$partial_fasta_tmp" > "$combined_fasta"

rm -f "$full_patterns" "$partial_patterns" "$filtered_headers"
rm -f "$full_fasta_tmp" "$partial_fasta_tmp"
rm -f "$full_contigs" "$partial_contigs"

exit 0
