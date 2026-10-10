#!/usr/bin/env python

"""Combine per-sample KneadData read counts into a cohort table and MultiQC sections.

Writes <prefix>_read_counts.tsv (long format, the analysis input), a MultiQC table,
and a MultiQC stacked bar graph. The bar graph is written only for a cohort small
enough to read and only under serial filtering, where per-database removals
partition the reads; under parallel filtering they overlap and must not be stacked.

See docs/workflows/qc_dna.md.
"""

import argparse
import csv
from pathlib import Path

STAGES = ("raw", "trimmed", "decontaminated", "final")


def read_counts(paths):
    """Load the per-sample TSVs into rows plus the per-sample run settings."""
    rows = []
    settings = {}

    for path in paths:
        sample_settings = {}
        with path.open() as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                key, _, value = line[1:].strip().partition("\t")
                sample_settings[key.strip()] = value.strip()

        with path.open() as handle:
            body = (line for line in handle if not line.startswith("#"))
            for row in csv.DictReader(body, delimiter="\t"):
                row["count"] = int(row["count"])
                rows.append(row)

        sample = sample_settings.get("sample")
        if sample:
            settings[sample] = sample_settings

    return rows, settings


def database_order(rows):
    """Databases in the order KneadData filtered against them.

    First-seen order in the log is the order of the chain under serial filtering,
    which is the order the stacked bar has to follow to read correctly.
    """
    order = []
    for row in rows:
        database = row["database"]
        if row["stage"] == "decontaminated" and database and database not in order:
            order.append(database)
    return order


def summarise(rows, databases, mate="pair1"):
    """One record per sample: raw, what each step removed, and what survived."""
    by_sample = {}
    for row in rows:
        if row["mate"] != mate:
            continue
        by_sample.setdefault(row["sample"], {})[(row["stage"], row["database"])] = row["count"]

    summary = {}
    for sample, counts in sorted(by_sample.items()):
        raw = counts.get(("raw", ""))
        trimmed = counts.get(("trimmed", ""))
        final = counts.get(("final", ""))
        if raw is None or trimmed is None or final is None:
            continue

        record = {"raw": raw, "trimmed": trimmed, "final": final, "lost_trimming": raw - trimmed}
        previous = trimmed
        for database in databases:
            surviving = counts.get(("decontaminated", database))
            if surviving is None:
                continue
            record[f"removed_{database}"] = previous - surviving
            record[f"surviving_{database}"] = surviving
            previous = surviving
        summary[sample] = record

    return summary


def summarise_parallel(rows, databases, mate="pair1"):
    """As `summarise`, but every database is measured against the trimmed reads."""
    summary = summarise(rows, databases, mate)
    by_sample = {}
    for row in rows:
        if row["mate"] == mate and row["stage"] == "decontaminated":
            by_sample.setdefault(row["sample"], {})[row["database"]] = row["count"]

    for sample, record in summary.items():
        for database in databases:
            surviving = by_sample.get(sample, {}).get(database)
            if surviving is not None:
                record[f"removed_{database}"] = record["trimmed"] - surviving
    return summary


def write_long_table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "stage", "database", "mate", "count"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_mqc_table(path, summary, databases, mode):
    columns = ["raw", "lost_trimming"] + [f"removed_{db}" for db in databases] + ["final"]
    header = [
        "# id: 'kneaddata_read_fate_table'",
        "# section_name: 'Trimming and decontamination'",
        f"# description: 'Read pairs at each stage, with one column per reference database. "
        f"KneadData filtered in {mode} mode.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'kneaddata_read_fate_table_plot'",
        "#     namespace: 'KneadData'",
        "#     col1_header: 'Sample'",
        "#     sort_rows: false",
    ]
    with path.open("w", newline="") as handle:
        handle.write("\n".join(header) + "\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample"] + columns)
        for sample, record in summary.items():
            writer.writerow([sample] + [record.get(column, 0) for column in columns])


def write_mqc_bargraph(path, summary, databases, mode):
    categories = ["final"] + [f"removed_{db}" for db in reversed(databases)] + ["lost_trimming"]
    labels = {
        "final": "Surviving",
        "lost_trimming": "Removed by trimming",
        **{f"removed_{db}": f"Matched {db}" for db in databases},
    }

    lines = [
        "id: 'kneaddata_read_fate'",
        "section_name: 'Read fate'",
        "description: >-",
        f"    Where every raw read pair went. KneadData filtered in {mode} mode, so the",
        "    parts sum to the raw read count and each database is credited only with",
        "    reads the databases before it had not already removed.",
        "plot_type: 'bargraph'",
        "pconfig:",
        "    id: 'kneaddata_read_fate_plot'",
        "    title: 'KneadData: read fate'",
        "    ylab: 'Read pairs'",
        "    cpswitch_counts_label: 'Read pairs'",
        "categories:",
    ]
    for category in categories:
        lines.append(f"    {category}:")
        lines.append(f"        name: '{labels[category]}'")
    lines.append("data:")
    for sample, record in summary.items():
        lines.append(f"    '{sample}':")
        for category in categories:
            lines.append(f"        {category}: {record.get(category, 0)}")

    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts", type=Path, nargs="+", help="Per-sample read-count TSVs")
    parser.add_argument("--prefix", default="kneaddata", help="Output file prefix")
    parser.add_argument(
        "--max-plot-samples",
        type=int,
        default=100,
        help="Draw the bar graph only up to this many samples (default: 100). "
        "Beyond it the table carries the same numbers and stays readable.",
    )
    args = parser.parse_args()

    rows, settings = read_counts(args.counts)
    if not rows:
        raise SystemExit("No read counts found in the supplied files.")

    modes = {s.get("filtering_mode", "unknown") for s in settings.values()}
    mode = modes.pop() if len(modes) == 1 else "mixed"

    databases = database_order(rows)
    summary = (
        summarise(rows, databases)
        if mode == "serial"
        else summarise_parallel(rows, databases)
    )

    write_long_table(Path(f"{args.prefix}_read_counts.tsv"), rows)
    write_mqc_table(Path(f"{args.prefix}_table_mqc.tsv"), summary, databases, mode)

    if mode != "serial":
        print(
            f"Filtering mode is '{mode}', so per-database removals overlap and cannot be "
            f"stacked. Wrote the table only; the per-database columns are still there."
        )
    elif len(summary) > args.max_plot_samples:
        print(
            f"{len(summary)} samples exceeds --max-plot-samples={args.max_plot_samples}; "
            f"wrote the table only. The full counts are in {args.prefix}_read_counts.tsv."
        )
    else:
        write_mqc_bargraph(Path(f"{args.prefix}_bargraph_mqc.yaml"), summary, databases, mode)


if __name__ == "__main__":
    main()
