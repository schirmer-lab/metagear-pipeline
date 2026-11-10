#!/usr/bin/env python3
"""Generate the mock metagenome used by the end-to-end test profile.

The bundled `test` profile checks that the pipeline installs and its plumbing
works; it runs no analysis workflow. This generates a small synthetic community
that the analysis workflows can process to completion in a few minutes, so an
end-to-end run can be exercised with no reference database and no cohort data.

Design constraints, in order:

1. Hermetic. No network, no external test-data repository, no reference
   database. Everything the test needs is produced here and committed.
2. Assemblable. MEGAHIT needs real overlap structure, so reads are simulated
   from defined source sequences at sufficient depth rather than drawn at
   random.
3. Gene-bearing. Prodigal must find coding sequences, so each source sequence
   carries planted open reading frames with a start codon, no in-frame stop,
   and a terminal stop. Without this, gene calling returns nothing and the
   clustering, annotation and abundance stages have no input.
4. Small. The FASTQs are committed, so the whole set stays a few MB.
5. Deterministic. A fixed seed, so a snapshot assertion is stable.

Usage:
    python tests/make_mock_metagenome.py --outdir assets/test_e2e
"""
import argparse
import gzip
import pathlib
import random

# Sense codons only, so a planted ORF never carries an in-frame stop.
CODONS = [
    "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "AAT", "AAC", "GAT", "GAC",
    "TGT", "TGC", "CAA", "CAG", "GAA", "GAG", "GGT", "GGC", "CAT", "CAC",
    "ATT", "ATC", "ATA", "TTA", "TTG", "CTT", "CTC", "AAA", "AAG", "ATG",
    "TTT", "TTC", "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "AGT", "AGC",
    "ACT", "ACC", "ACA", "ACG", "TGG", "TAT", "TAC", "GTT", "GTC", "GTA",
]
STOPS = ("TAA", "TAG", "TGA")


def orf(rng, n_codons):
    body = "".join(rng.choice(CODONS) for _ in range(n_codons))
    return "ATG" + body + rng.choice(STOPS)


def intergenic(rng, length):
    """AT-rich spacer, so Prodigal does not run adjacent genes together."""
    return "".join(rng.choices("ATATATGC", k=length))


def make_genome(rng, n_genes, gene_codons=(80, 320), spacer=(60, 200)):
    parts = []
    for _ in range(n_genes):
        parts.append(orf(rng, rng.randint(*gene_codons)))
        parts.append(intergenic(rng, rng.randint(*spacer)))
    return "".join(parts)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def mutate(rng, seq, err):
    if err <= 0:
        return seq
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < err:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def simulate(rng, genome, depth, read_len, frag_mean, frag_sd, err):
    n_pairs = max(1, int(len(genome) * depth / (2 * read_len)))
    r1, r2 = [], []
    for _ in range(n_pairs):
        frag = max(read_len + 20, int(rng.gauss(frag_mean, frag_sd)))
        frag = min(frag, len(genome))
        start = rng.randint(0, len(genome) - frag)
        piece = genome[start:start + frag]
        r1.append(mutate(rng, piece[:read_len], err))
        r2.append(mutate(rng, revcomp(piece[-read_len:]), err))
    return r1, r2


def write_fastq(path, reads, prefix):
    with gzip.open(path, "wt", compresslevel=9) as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@{prefix}_{i}/1\n{seq}\n+\n{'I' * len(seq)}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--samples", type=int, default=3)
    ap.add_argument("--genomes", type=int, default=3)
    ap.add_argument("--genes-per-genome", type=int, default=12)
    ap.add_argument("--depth", type=float, default=22.0)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=20261019)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    outdir = args.outdir
    reads_dir = outdir / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)

    # One shared community across samples, at rotating relative abundance, so
    # cross-sample gene clustering and abundance estimation have real signal.
    genomes = {f"mock_genome_{i + 1}": make_genome(rng, args.genes_per_genome)
               for i in range(args.genomes)}
    with (outdir / "mock_community.fasta").open("w") as fh:
        for name, seq in genomes.items():
            fh.write(f">{name} planted_orfs={args.genes_per_genome} len={len(seq)}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    rows = ["sample,biome,fastq_1,fastq_2"]
    for s in range(args.samples):
        sid = f"MOCK{s + 1}"
        r1_all, r2_all = [], []
        for gi, (_, seq) in enumerate(genomes.items()):
            weight = [1.0, 0.6, 0.3][(gi + s) % 3]
            r1, r2 = simulate(rng, seq, args.depth * weight, args.read_len,
                              frag_mean=320, frag_sd=40, err=0.002)
            r1_all += r1
            r2_all += r2
        order = list(range(len(r1_all)))
        rng.shuffle(order)
        r1_all = [r1_all[i] for i in order]
        r2_all = [r2_all[i] for i in order]
        write_fastq(reads_dir / f"{sid}_R1.fastq.gz", r1_all, f"{sid}_read")
        write_fastq(reads_dir / f"{sid}_R2.fastq.gz", r2_all, f"{sid}_read")
        # __TEST_E2E_DIR__ is substituted by conf/test_e2e.config at parse time,
        # so the committed samplesheet carries no machine-specific path.
        # biome must be a SemiBin2 pretrained environment or "global" (schema-validated)
        rows.append(f"{sid},global,__TEST_E2E_DIR__/reads/{sid}_R1.fastq.gz,"
                    f"__TEST_E2E_DIR__/reads/{sid}_R2.fastq.gz")
        print(f"{sid}: {len(r1_all)} read pairs")

    (outdir / "samplesheet.csv").write_text("\n".join(rows) + "\n")

    # Stand-ins that let the profile skip protein annotation; see docs/end_to_end_test.md
    (outdir / "stub_protein_annotations.tsv").write_text(
        "# Stub standing in for PROTEIN_ANNOTATION output; see docs/end_to_end_test.md\n"
        "protein_id\tannotation\n")
    amr = outdir / "stub_amrfinder_db"
    amr.mkdir(exist_ok=True)
    (amr / "README.txt").write_text(
        "Placeholder. GENES_INIT requires params.amrfinder_db to exist; nothing here is read.\n")

    total = sum(p.stat().st_size for p in reads_dir.glob("*.fastq.gz"))
    print(f"wrote {outdir}/samplesheet.csv; reads total {total / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
