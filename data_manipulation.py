#!/usr/bin/env python3

"""
Filter a protein2ipr mapping to include only Swiss-Prot accessions.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import IO, Iterable, Iterator, Optional, Set

UNIPROT_ACC_RE = re.compile(
    # Based on UniProt accession format; reused to validate parsed IDs.
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a protein2ipr mapping file to include only Swiss-Prot proteins. "
            "Supports plain text or gzip-compressed inputs."
        )
    )
    parser.add_argument(
        "--protein2ipr",
        required=True,
        type=Path,
        help="Path to the protein2ipr.dat (optionally .gz) mapping file.",
    )
    parser.add_argument(
        "--swiss-prot-fasta",
        required=True,
        type=Path,
        help="Path to the Swiss-Prot FASTA file containing the target accessions.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output path for the filtered mapping (plain text or .gz).",
    )
    return parser.parse_args()


def open_text_file(path: Path, mode: str) -> IO[str]:
    if "t" not in mode:
        raise ValueError("open_text_file should be used with text modes only.")
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="ignore")
    return open(path, mode, encoding="utf-8", errors="ignore")


def extract_accession_from_header(header: str) -> Optional[str]:
    primary_token = header.split()[0] if header else ""
    if primary_token.startswith("sp|"):
        parts = primary_token.split("|")
        if len(parts) >= 2 and parts[1]:
            accession = parts[1]
            if UNIPROT_ACC_RE.match(accession):
                return accession
    for token in re.split(r"[| ,;]+", header):
        if UNIPROT_ACC_RE.match(token):
            return token
    return None


def load_swiss_prot_accessions(fasta_path: Path) -> Set[str]:
    accessions: Set[str] = set()
    with open_text_file(fasta_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                accession = extract_accession_from_header(header)
                if accession:
                    accessions.add(accession)
    return accessions


def read_protein2ipr(mapping_path: Path) -> Iterator[str]:
    with open_text_file(mapping_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            yield line


def filter_protein2ipr(
    lines: Iterable[str], swiss_accessions: Set[str]
) -> Iterator[str]:
    for line in lines:
        accession = line.split(None, 1)[0]
        if accession in swiss_accessions:
            yield line


def write_filtered_database(lines: Iterable[str], output_path: Path) -> int:
    count = 0
    with open_text_file(output_path, "wt") as handle:
        for line in lines:
            handle.write(line)
            count += 1
    return count


def main() -> None:
    args = parse_args()
    swiss_accessions = load_swiss_prot_accessions(args.swiss_prot_fasta)
    if not swiss_accessions:
        print(
            "No Swiss-Prot accessions were detected in the FASTA file.",
            file=sys.stderr,
        )
        sys.exit(1)

    input_lines = read_protein2ipr(args.protein2ipr)
    filtered_lines = filter_protein2ipr(input_lines, swiss_accessions)
    written = write_filtered_database(filtered_lines, args.output)
    print(f"Wrote {written} Swiss-Prot records to {args.output}")


if __name__ == "__main__":
    main()
