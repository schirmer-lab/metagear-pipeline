# Per-contig taxonomy enrichment

After `integrated_classification`, each sample's `classification/per_contig/<sample>.contigs.tsv`
attributes contigs to Binette MAGs but does **not** carry species-level lineage —
`integrated_classification` only knows a contig is in a MAG, not what species that MAG is.

`cohort_dereplication` produces all the pieces needed to fill in the lineage:

- **`cohort_dereplication/drep/data_tables/Cdb.csv`** — per-genome cluster table (each MAG → secondary cluster)
- **`cohort_dereplication/drep/data_tables/Wdb.csv`** — cluster winners (cluster → representative genome)
- **`cohort_dereplication/gtdbtk/gtdbtk.{bac120,ar53}.summary.tsv`** — GTDB-Tk classification of each cluster representative
- **`assemblies/bins/<sample>/<sample>.contig_to_bin.tsv`** — per-sample contig → bin mapping (from Binette)
- **`classification/per_contig/<sample>.contigs.tsv`** — per-sample per-contig classification (from MERGE_CONTIG_CLASSIFICATION)

The join is intentionally **not** in the pipeline — it's a few-line pandas operation
that runs in seconds, and the same building blocks support several other join shapes
(contig × MAG lineage, contig × MAG abundance, MAG × abundance × lineage, etc.).

## Pandas recipe (per sample)

```python
import pandas as pd
from pathlib import Path

# ─── Paths ─────────────────────────────────────────────────────────────────
RESULTS = Path("/path/to/results")
SAMPLE  = "SAMPLE-0"

per_contig_tsv = RESULTS / f"classification/per_contig/{SAMPLE}.contigs.tsv"
c2b_tsv        = RESULTS / f"assemblies/bins/{SAMPLE}/{SAMPLE}.contig_to_bin.tsv"
cdb_csv        = RESULTS / "cohort_dereplication/drep/data_tables/Cdb.csv"
wdb_csv        = RESULTS / "cohort_dereplication/drep/data_tables/Wdb.csv"
gtdb_summaries = sorted((RESULTS / "cohort_dereplication/gtdbtk").glob("gtdbtk.*.summary.tsv"))

# ─── Load ──────────────────────────────────────────────────────────────────
per_contig = pd.read_csv(per_contig_tsv, sep="\t")
c2b        = pd.read_csv(c2b_tsv,        sep="\t", names=["contig_id", "bin_name"])
cdb        = pd.read_csv(cdb_csv)
wdb        = pd.read_csv(wdb_csv)
gtdb       = pd.concat([pd.read_csv(f, sep="\t") for f in gtdb_summaries], ignore_index=True)

# ─── Join chain ────────────────────────────────────────────────────────────
# bin_name (Binette) -> dRep genome filename (Cdb.genome) -> secondary cluster ->
# cluster winner (Wdb.genome) -> GTDB lineage (gtdb.classification)
c2b["bin_file"] = c2b["bin_name"] + ".fa"
out = (
    per_contig
      .merge(c2b, on="contig_id", how="left")
      .merge(cdb[["genome", "secondary_cluster"]],
             left_on="bin_file", right_on="genome", how="left")
      .merge(wdb[["cluster", "genome"]].rename(columns={"genome": "winner_file"}),
             left_on="secondary_cluster", right_on="cluster", how="left")
      .assign(winner_id=lambda d: d["winner_file"].str.replace(r"\.fa$", "", regex=True))
      .merge(gtdb[["user_genome", "classification"]].rename(columns={"classification": "gtdb_lineage"}),
             left_on="winner_id", right_on="user_genome", how="left")
)

# ─── Result columns added ──────────────────────────────────────────────────
# - secondary_cluster: dRep cluster id (e.g. "1_2")
# - winner_id:         representative MAG name
# - gtdb_lineage:      full GTDB-Tk lineage (d__Bacteria;p__...;s__...)

out.to_csv(f"{SAMPLE}.contigs.enriched.tsv", sep="\t", index=False)
```

## Cohort-wide version

For all 20 samples at once:

```python
import pandas as pd
from pathlib import Path

RESULTS = Path("/path/to/results")
cdb     = pd.read_csv(RESULTS / "cohort_dereplication/drep/data_tables/Cdb.csv")
wdb     = pd.read_csv(RESULTS / "cohort_dereplication/drep/data_tables/Wdb.csv")
gtdb    = pd.concat([pd.read_csv(f, sep="\t")
                     for f in (RESULTS / "cohort_dereplication/gtdbtk").glob("gtdbtk.*.summary.tsv")],
                    ignore_index=True)

# Pre-build a single bin_file -> gtdb_lineage map for the whole cohort
bin_to_lineage = (
    cdb[["genome", "secondary_cluster"]]
      .merge(wdb[["cluster", "genome"]].rename(columns={"genome": "winner_file"}),
             left_on="secondary_cluster", right_on="cluster")
      .assign(winner_id=lambda d: d["winner_file"].str.replace(r"\.fa$", "", regex=True))
      .merge(gtdb[["user_genome", "classification"]].rename(columns={"classification": "gtdb_lineage"}),
             left_on="winner_id", right_on="user_genome", how="left")
      [["genome", "secondary_cluster", "gtdb_lineage"]]
)

# Apply to each sample's per-contig TSV
for per_contig_tsv in (RESULTS / "classification/per_contig").glob("*.contigs.tsv"):
    sample = per_contig_tsv.stem.replace(".contigs", "")
    c2b    = pd.read_csv(RESULTS / f"assemblies/bins/{sample}/{sample}.contig_to_bin.tsv",
                         sep="\t", names=["contig_id", "bin_name"])
    c2b["bin_file"] = c2b["bin_name"] + ".fa"

    enriched = (
        pd.read_csv(per_contig_tsv, sep="\t")
          .merge(c2b, on="contig_id", how="left")
          .merge(bin_to_lineage, left_on="bin_file", right_on="genome", how="left")
    )
    enriched.to_csv(f"{sample}.contigs.enriched.tsv", sep="\t", index=False)
```

## DuckDB version (very large cohorts)

For thousands of samples / millions of contigs, where loading everything in pandas
is awkward:

```sql
CREATE TABLE cdb  AS SELECT * FROM read_csv_auto('results/cohort_dereplication/drep/data_tables/Cdb.csv');
CREATE TABLE wdb  AS SELECT * FROM read_csv_auto('results/cohort_dereplication/drep/data_tables/Wdb.csv');
CREATE TABLE gtdb AS SELECT user_genome, classification
                    FROM read_csv_auto('results/cohort_dereplication/gtdbtk/gtdbtk.*.summary.tsv',
                                       union_by_name=true);

-- bin file -> gtdb lineage lookup (cohort-wide)
CREATE TABLE bin_to_lineage AS
SELECT cdb.genome     AS bin_file,
       cdb.secondary_cluster,
       gtdb.classification AS gtdb_lineage
FROM   cdb
       JOIN wdb  ON wdb.cluster      = cdb.secondary_cluster
       LEFT JOIN gtdb ON gtdb.user_genome = REGEXP_REPLACE(wdb.genome, '\.fa$', '');

-- Per-sample enrichment (loop over samples or use a wildcard read)
COPY (
    WITH c2b AS (
        SELECT column0 AS contig_id, column1 || '.fa' AS bin_file
        FROM read_csv_auto('results/assemblies/bins/SAMPLE-0/SAMPLE-0.contig_to_bin.tsv',
                           header=false, delim='\t')
    )
    SELECT pc.*, c2b.bin_file, btl.secondary_cluster, btl.gtdb_lineage
    FROM read_csv_auto('results/classification/per_contig/SAMPLE-0.contigs.tsv') pc
         LEFT JOIN c2b           ON c2b.contig_id = pc.contig_id
         LEFT JOIN bin_to_lineage btl ON btl.bin_file = c2b.bin_file
) TO 'SAMPLE-0.contigs.enriched.tsv' (HEADER, DELIMITER '\t');
```

## What you get

Per row of `<sample>.contigs.tsv`, the enrichment adds:

| Column | Source | Filled when |
|---|---|---|
| `bin_file` | `<sample>.contig_to_bin.tsv` (Binette) | contig is in a MAG |
| `secondary_cluster` | `Cdb.csv` (dRep) | bin is in a species cluster |
| `winner_id` (optional) | `Wdb.csv` (dRep) | cluster has a representative |
| `gtdb_lineage` | GTDB-Tk summary | cluster winner was successfully placed |

For unbinned chromosome contigs, all four columns stay null — those contigs already
have whatever signal MERGE_CONTIG_CLASSIFICATION provided (Tiara label, etc.).

## See also

- `docs/workflows/integrated_classification.md` — per-contig TSV schema
- `docs/workflows/cohort_dereplication.md` — cohort outputs (dRep tables, GTDB-Tk
  summaries, MAG catalog)
