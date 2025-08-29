process VIRSORTER2_SETUP {
  label 'process_single'

  // Use either container or conda; container takes precedence when both are set by your config
  conda "bioconda::virsorter=2.2.4"   // version pin (build tag not supported here)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/virsorter:2.2.4--pyhdfd78af_0' :
      'quay.io/biocontainers/virsorter:2.2.4--pyhdfd78af_0' }"

  output:
  path "virsorter2_db/", emit: virsorter2_db
  path "versions.yml" ,  emit: versions

  script:
  """
  set -euo pipefail

  # Write a wget→curl shim into the task work dir (handles --opt=VAL and --opt VAL)
  cat > wget <<'EOF'
  #!/bin/sh
  OUTFILE=""
  URL=""
  RETRIES=2
  WAITRETRY=60
  TIMEOUT=60

  while [ \$# -gt 0 ]; do
    case "\$1" in
      -O) OUTFILE="\$2"; shift 2 ;;
      -O=*) OUTFILE="\${1#*=}"; shift ;;
      -q|-nv|--quiet|--no-verbose) shift ;;
      --tries) RETRIES="\$2"; shift 2 ;;
      --tries=*) RETRIES="\${1#*=}"; shift ;;
      --retry-connrefused) shift ;;   # curl will use --retry-all-errors
      --waitretry) WAITRETRY="\$2"; shift 2 ;;
      --waitretry=*) WAITRETRY="\${1#*=}"; shift ;;
      --timeout) TIMEOUT="\$2"; shift 2 ;;
      --timeout=*) TIMEOUT="\${1#*=}"; shift ;;
      http://*|https://*)
        [ -z "\$URL" ] && URL="\$1"; shift ;;
      *)
        # Ignore other wget-only flags to avoid leaking to curl
        shift ;;
    esac
  done

  [ -z "\$OUTFILE" ] && OUTFILE="/dev/stdout"
  [ -z "\$URL" ] && { echo "wget-shim: no URL" >&2; exit 2; }

  TOTAL_MAX=\$(( (RETRIES + 1) * TIMEOUT * 2 ))

  echo "[wget-shim] curl -L -f -s -S --retry \$RETRIES --retry-all-errors --retry-delay \$WAITRETRY --connect-timeout \$TIMEOUT --max-time \$TOTAL_MAX -o \$OUTFILE \$URL" >&2

  exec curl -L -f -s -S \
    --retry "\$RETRIES" --retry-all-errors --retry-delay "\$WAITRETRY" \
    --connect-timeout "\$TIMEOUT" --max-time "\$TOTAL_MAX" \
    -o "\$OUTFILE" "\$URL"
  EOF
  chmod +x wget

  # Prove the shim is first on PATH for THIS command invocation
  echo "PATH before: \$PATH"
  PATH="\$PWD:\$PATH"
  echo "PATH after : \$PATH"
  type -a wget || true
  ./wget -O /dev/null https://example.com || true

  # Run setup with shim in effect
  PATH="\$PWD:\$PATH" virsorter setup -d virsorter2_db -j ${task.cpus}

  # Record versions
  cat <<-END_YAML > versions.yml
  "${task.process}":
    VirSorter2: \$(grep -m1 'VirSorter' .command.log | awk '{print \$NF}')
  END_YAML
  """
}