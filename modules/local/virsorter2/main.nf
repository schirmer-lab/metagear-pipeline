process VIRSORTER2 {
  tag "$meta.id"

  conda "bioconda::virsorter==2.2.4--pyhdfd78af_1"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/virsorter:2.2.4--pyhdfd78af_0' :
      'biocontainers/virsorter:2.2.4--pyhdfd78af_0' }"

  input:
  tuple val(meta), path(contigs)
  path db

  output:
  tuple val(meta), path("*/*.final-viral-combined.fa"), emit: vs2_virus
  tuple val(meta), path("*/*.final-viral-score.tsv"),   emit: vs2_score
  tuple val(meta), path("*/for-dramv/*.final-viral-combined-for-dramv.fa"), optional: true, emit: vs2_4dram_virus
  tuple val(meta), path("*/for-dramv/*.viral-affi-contigs-for-dramv.tab"),  optional: true, emit: vs2_4dra_affi
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  // User-tunable extras
  def prefix  = task.ext.prefix ?: "${meta.id}"
  def args    = task.ext.args   ?: ''   // e.g., --keep-original-seq ...
  def args2   = task.ext.args2  ?: ''   // optional suffix for workdir name
  def vs2_dir = "${prefix}${args2}"

  """
  set -euo pipefail

  # ---------- Local caches (per task, node-local) ----------
  mkdir -p "\$PWD/.conda_pkgs" "\$PWD/_tmp" "\$PWD/.vs2_conda"
  export CONDA_PKGS_DIRS="\$PWD/.conda_pkgs"
  export TMPDIR="\$PWD/_tmp"
  export SNAKEMAKE_CONDA_PREFIX="\$PWD/.vs2_conda"
  umask 0002

  # ---------- Input handling ----------
  INPUT="${contigs}"
  if [[ "${contigs}" == *.gz ]]; then
      gunzip -c "${contigs}" >| "${prefix}.decompressed.fasta"
      INPUT="${prefix}.decompressed.fasta"
  fi

  # ---------- Run VirSorter2 (use real DB; override conda prefix) ----------
  virsorter run ${args} -i "\$INPUT" -w "${vs2_dir}" -j ${task.cpus} -d "${db}" all \
      --conda-prefix "\$PWD/.vs2_conda" \
      --conda-frontend mamba \
      --latency-wait 600 \
      --rerun-incomplete \
      --use-conda

  # ---------- Normalize output names ----------
  mv "${vs2_dir}/final-viral-combined.fa"            "${vs2_dir}/${prefix}.final-viral-combined.fa"
  mv "${vs2_dir}/final-viral-score.tsv"              "${vs2_dir}/${prefix}.final-viral-score.tsv"
  mv "${vs2_dir}/for-dramv/final-viral-combined-for-dramv.fa" "${vs2_dir}/for-dramv/${prefix}.final-viral-combined-for-dramv.fa" || true
  mv "${vs2_dir}/for-dramv/viral-affi-contigs-for-dramv.tab"  "${vs2_dir}/for-dramv/${prefix}.viral-affi-contigs-for-dramv.tab"  || true

  # Optional cleanup
  [[ -f "${prefix}.decompressed.fasta" ]] && rm -f "${prefix}.decompressed.fasta" || true

  # ---------- Provenance (robust to set -e) ----------
  VS2_VER="\$(virsorter --version 2>/dev/null | awk '{print \$NF}' || echo unknown)"
  cat <<-END_VERSIONS > versions.yml
  "METAGEAR:VIRUS:VIRAL_ANNOTATION:VIRSORTER2${args2}":
    VirSorter2: \${VS2_VER}
  END_VERSIONS
  """
}
