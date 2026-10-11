#!/usr/bin/env python

"""Turn a KneadData log into a per-sample table of read counts.

KneadData counts reads at every stage, and once per reference database when it
decontaminates, but writes those counts only to its log. This reads them back out.

Stages and the serial-vs-parallel caveat: docs/workflows/qc_dna.md.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

READ_COUNT = re.compile(
    r"READ COUNT:\s+"
    r"(?P<stage>raw|trimmed|decontaminated|final)\s+"
    r"(?:(?P<database>.+?)\s+)?"
    r"(?P<mate>pair1|pair2|orphan1|orphan2|single)\s+:"
    r".*:\s*(?P<count>[0-9]+(?:\.[0-9]+)?)\s*$"
)
ARGUMENT = re.compile(r"^(?P<key>[a-z_0-9]+)\s*=\s*(?P<value>.*)$")
VERSION = re.compile(r"Running kneaddata v(?P<version>\S+)")

FIELDNAMES = ["sample", "stage", "database", "mate", "count"]


def parse_log(handle):
    """Read one KneadData log into (counts, run settings)."""
    counts = []
    settings = {}
    in_arguments = False

    for line in handle:
        version = VERSION.search(line)
        if version:
            settings["kneaddata_version"] = version.group("version")

        if "Running with the following arguments" in line:
            in_arguments = True
            continue

        count = READ_COUNT.search(line)
        if count:
            in_arguments = False
            counts.append(
                {
                    "stage": count.group("stage"),
                    "database": (count.group("database") or "").strip(),
                    "mate": count.group("mate"),
                    # KneadData writes these as floats ("91515547.0").
                    "count": int(float(count.group("count"))),
                }
            )
            continue

        if in_arguments:
            argument = ARGUMENT.match(line.strip())
            if argument:
                settings[argument.group("key")] = argument.group("value").strip()
            elif line.strip() and " - " in line:
                in_arguments = False

    return counts, settings


def filtering_mode(settings):
    """serial or parallel, as the log reports it.

    KneadData coerces --serial on before it logs its arguments (0.12 forces it
    whenever --processes is 1), so the logged value is the effective one and no
    inference from the version or the process count is needed.
    """
    serial = settings.get("serial")
    if serial is None:
        return "unknown"
    return "serial" if serial.lower() == "true" else "parallel"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, help="KneadData log file")
    parser.add_argument("output", type=Path, help="TSV to write")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    args = parser.parse_args()

    with args.log.open() as handle:
        counts, settings = parse_log(handle)

    if not counts:
        sys.exit(f"No READ COUNT lines found in {args.log}. KneadData did not finish.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        handle.write(f"# sample\t{args.sample}\n")
        handle.write(f"# kneaddata_version\t{settings.get('kneaddata_version', 'unknown')}\n")
        handle.write(f"# filtering_mode\t{filtering_mode(settings)}\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in counts:
            writer.writerow({"sample": args.sample, **row})


if __name__ == "__main__":
    main()
