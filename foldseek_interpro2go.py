import time
import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import tempfile
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Set
from GOTool import GeneOntology


# Utils for parsing InterPro2GO and protein2ipr.dat
UNIPROT_ACC_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])$"
)
IPR_RE = re.compile(r"(IPR\d{5,})")

FASTA_HEADER_RE = re.compile(r"^>(.*)")
NON_AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")


def extract_uniprot_from_id(identifier: str) -> Optional[str]:
    """
    Extract a UniProt accession from Foldseek target identifier
    (AlphaFold/UniProt styles supported).
    """
    m = re.search(r"AF-([A-Z0-9]+)-F\d", identifier)  # AlphaFold style
    if m:
        return m.group(1)
    m = re.search(r"(?:sp|tr)\|([A-Z0-9]+)\|", identifier)  # UniProt FASTA style
    if m:
        return m.group(1)
    for tok in re.split(
        r"\W+", identifier
    ):  # fallback: any token that looks like an accession
        if UNIPROT_ACC_RE.match(tok):
            return tok
    return None


def open_maybe_gzip(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def parse_interpro2go(path: Optional[str]) -> Dict[str, Set[str]]:
    if not path:
        return {}
    ipr2go: Dict[str, Set[str]] = defaultdict(set)
    with open_maybe_gzip(path) as f:
        for line in f:
            if not line.strip() or line.startswith("!"):
                continue
            m_ipr = re.search(r"InterPro:(IPR\d{5,})", line)
            if not m_ipr:
                continue
            ipr = m_ipr.group(1)
            for m_go in re.finditer(r"(GO:\d{7})", line):
                ipr2go[ipr].add(m_go.group(1))
    return ipr2go


def stream_protein2ipr_for_wanted(
    mapping_path: str, wanted: Set[str]
) -> Dict[str, Set[str]]:
    """
    Stream protein2ipr.dat[.gz] and
    return mappings only for accessions in `wanted`.
    """
    acc2ipr: Dict[str, Set[str]] = defaultdict(set)
    if not wanted:
        return {}
    with open_maybe_gzip(mapping_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if not parts:
                continue
            acc = parts[0]
            if acc in wanted:
                for ipr in IPR_RE.findall(line):
                    acc2ipr[acc].add(ipr)
    return acc2ipr


@dataclass
class FastaSeq:
    header: str
    seq: str
    name: str  # sanitized name for file basenames


def read_fasta(path: Path) -> List[FastaSeq]:
    """
    Minimal FASTA reader returning a list of sequences with sanitized names.
    """
    records: List[FastaSeq] = []
    header = None
    seq_chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = re.sub(r"\s+", "", "".join(seq_chunks)).upper()
                    seq = NON_AA_RE.sub("", seq)
                    name = sanitize_name(header)
                    records.append(FastaSeq(header=header, seq=seq, name=name))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            seq = re.sub(r"\s+", "", "".join(seq_chunks)).upper()
            seq = NON_AA_RE.sub("", seq)
            name = sanitize_name(header)
            records.append(FastaSeq(header=header, seq=seq, name=name))
    return records


def write_single_fasta(out_path: Path, header: str, seq: str) -> None:
    out_path.write_text(f">{header}\n{wrap_seq(seq)}\n", encoding="utf-8")


def wrap_seq(s: str, width: int = 60) -> str:
    return "\n".join(s[i : i + width] for i in range(0, len(s), width))  # noqa


def sanitize_name(header: str) -> str:
    # Use first token before spaces and strip non-word characters
    base = header.split()[0]
    base = re.sub(r"[^\w\-\.]+", "_", base)
    return base[:80]  # keep it manageable


# ------------------------------- Foldseek -----------------------------------


@dataclass
class Hit:
    query: str
    target: str
    evalue: float
    bits: float
    fident: float
    alnlen: int
    qstart: int
    qend: int
    tstart: int
    tend: int
    qlen: int
    tlen: int
    qtmscore: Optional[float] = None
    ttmscore: Optional[float] = None
    alntmscore: Optional[float] = None
    lddt: Optional[float] = None

    @property
    def qcov(self) -> float:
        L = max(0, self.qend - self.qstart + 1)
        return L / self.qlen if self.qlen else 0.0

    @property
    def tcov(self) -> float:
        L = max(0, self.tend - self.tstart + 1)
        return L / self.tlen if self.tlen else 0.0

    @property
    def avg_tm(self) -> Optional[float]:
        if self.qtmscore is None or self.ttmscore is None:
            return None
        return 0.5 * (self.qtmscore + self.ttmscore)


def run_foldseek_easy_search(
    query_path: Path,
    target_db: str,
    tmpdir: Path,
    alignment_type: int,
    max_hits: int,
    prostt5_model: Optional[str] = None,
    gpu: Optional[str] = None,
) -> List[Hit]:
    """
    Run Foldseek easy-search and parse TSV results into Hit objects.

    Modes:
      - Structure query (PDB/mmCIF): default;
        can use --alignment-type 1 (TM-align re-score/rank).
      - FASTA query with ProstT5 (no coordinates): set prostt5_model to path;
        TM-align/TM-scores unavailable.
    """
    foldseek = shutil.which("foldseek")
    if not foldseek:
        raise RuntimeError("Foldseek not found in PATH.")

    out_tsv = tmpdir / f"{query_path.stem}.foldseek.tsv"
    FMT_FASTA = ",".join(
        [
            "query",
            "target",
            "evalue",
            "bits",
            "fident",
            "alnlen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "qlen",
            "tlen",
        ]
    )
    FMT_STRUCT = ",".join(
        [
            "query",
            "target",
            "evalue",
            "bits",
            "fident",
            "alnlen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "qlen",
            "tlen",
            "qtmscore",
            "ttmscore",
            "alntmscore",
            "lddt",
        ]
    )

    fmt = FMT_FASTA if prostt5_model else FMT_STRUCT

    cmd = [
        foldseek,
        "easy-search",
        str(query_path),
        target_db,
        str(out_tsv),
        str(tmpdir),
        "--format-output",
        fmt,
        "--max-seqs",
        str(max_hits),
    ]

    # FASTA→ProstT5 mode
    if prostt5_model:
        cmd += ["--prostt5-model", prostt5_model]
        # TM-align needs query Cα coordinates; disable it in ProstT5 mode.
        # (Keep alignment_type at 0; qtmscore/ttmscore/lddt will be empty.)
        alignment_type = 0
        cmd += ["--alignment-type", str(alignment_type)]

    # Structure→Structure mode: allow TM-align re-score if requested
    if alignment_type in (0, 1) and not prostt5_model:
        cmd += ["--alignment-type", str(alignment_type)]

    # GPU option if any
    if gpu:
        cmd += ["--gpu", gpu]

    print("[DEBUG] Running Foldseek command:", " ".join(cmd))
    cp = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )  # noqa
    if cp.returncode != 0:
        raise RuntimeError(
            f"Foldseek failed for {query_path.name}:\n{cp.stderr}"
        )  # noqa

    hits: List[Hit] = []
    with open(out_tsv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            def to_float(x: str) -> Optional[float]:
                try:
                    return float(x)
                except Exception:
                    return None

            hits.append(
                Hit(
                    query=parts[0],
                    target=parts[1],
                    evalue=float(parts[2]),
                    bits=float(parts[3]),
                    fident=float(parts[4]),
                    alnlen=int(float(parts[5])),
                    qstart=int(float(parts[6])),
                    qend=int(float(parts[7])),
                    tstart=int(float(parts[8])),
                    tend=int(float(parts[9])),
                    qlen=int(float(parts[10])),
                    tlen=int(float(parts[11])),
                    qtmscore=to_float(parts[12]) if len(parts) > 12 else None,
                    ttmscore=to_float(parts[13]) if len(parts) > 13 else None,
                    alntmscore=to_float(parts[14]) if len(parts) > 14 else None,  # noqa
                    lddt=to_float(parts[15]) if len(parts) > 15 else None,
                )
            )
    return hits


def select_hits(
    hits: List[Hit],
    evalue_max: float,
    min_qcov: float,
    min_tcov: float,
    min_avg_tm: Optional[float],
    topk: int,
    tm_available: bool,
) -> List[Hit]:
    sel: List[Hit] = []
    seen_targets: Set[str] = set()
    for h in hits:
        if h.evalue > evalue_max:
            continue
        if h.qcov < min_qcov or h.tcov < min_tcov:
            continue
        if tm_available and min_avg_tm is not None:
            if h.avg_tm is None or h.avg_tm < min_avg_tm:
                continue
        if h.target in seen_targets:
            continue
        sel.append(h)
        seen_targets.add(h.target)
        if len(sel) >= topk:
            break
    return sel


def group_hits_by_query(hits):
    byq = defaultdict(list)
    for h in hits:
        byq[h.query].append(h)
    return byq


# ----------------------------- Structure Handling ---------------------------
UNIPROT_ACC_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])$"
)


def extract_uniprot_from_header(header: str) -> Optional[str]:
    """
    Extract a UniProt accession from a FASTA header such as:
      sp_A0A0F6B8A2_BCSE_SALT1
      tr_A0A123ABC4_SOME_PROT
      >sp|Q3ZBH0|...
    Falls back to any UniProt-like token if the "sp_/tr_" form is absent.
    """
    # sp_/tr_ prefix like sp_A0A0F6B8A2_...
    m = re.search(r"(?:^|[>\s_])(?:sp|tr)_([A-Z0-9]+)", header)
    if m:
        return m.group(1)
    # UniProt FASTA form sp|ACC|...
    m = re.search(r"(?:sp|tr)\|([A-Z0-9]+)\|", header)
    if m:
        return m.group(1)
    # Fallback: any UniProt-like token
    for tok in re.split(r"[^\w]+", header):
        if UNIPROT_ACC_RE.match(tok):
            return tok
    return None


def normalise_query_identifier(query_id: str) -> str:
    """
    Attempt to map Foldseek query identifiers back to the protein identifiers
    used in downstream S2F processing (typically UniProt accessions).
    """
    if not query_id:
        return query_id

    cand = extract_uniprot_from_id(query_id)
    if cand:
        return cand

    cand = extract_uniprot_from_header(query_id)
    if cand:
        return cand

    # Tokens separated by common delimiters may still contain an accession.
    for tok in re.split(r"[^\w\-\|]+", query_id):
        cand = extract_uniprot_from_id(tok)
        if cand:
            return cand
        cand = extract_uniprot_from_header(tok)
        if cand:
            return cand

    return query_id


def ensure_structure_for_fasta(
    fasta_path: Path, mode: str, structures_dir: Path, search_recursive: bool = False
) -> List[Path]:
    """
    Ensure a structure file exists for each sequence in a FASTA.

    Supports:
      - Exact basename matches: {name}.pdb/.cif/.pdb.gz/.cif.gz
      - AlphaFold names: AF-<ACC>-F*-model*.{pdb|cif}[.gz]

    If search_recursive=True, use rglob under structures_dir.
    """
    seqs = read_fasta(fasta_path)
    out_paths: List[Path] = []
    globber = structures_dir.rglob if search_recursive else structures_dir.glob

    if mode == "existing":
        for rec in seqs:
            # 1) Try exact basename (four suffixes)
            candidates = [
                structures_dir / f"{rec.name}.pdb",
                structures_dir / f"{rec.name}.cif",
                structures_dir / f"{rec.name}.pdb.gz",
                structures_dir / f"{rec.name}.cif.gz",
            ]
            hit = next((p for p in candidates if p.exists()), None)

            # 2) If not found, derive UniProt accession from header and try AlphaFold pattern
            if hit is None:
                acc = extract_uniprot_from_header(rec.header)
                if acc:
                    # AlphaFold naming patterns (pdb/cif; gz/uncompressed; any model version)
                    patterns = [
                        f"AF-{acc}-F*-model*.pdb",
                        f"AF-{acc}-F*-model*.cif",
                        f"AF-{acc}-F*-model*.pdb.gz",
                        f"AF-{acc}-F*-model*.cif.gz",
                    ]
                    for pat in patterns:
                        found = sorted(globber(pat))
                        if found:
                            hit = found[0]
                            break

            if hit is None:
                raise FileNotFoundError(
                    f"No structure found for {rec.name}. "
                    f"Searched exact basenames and AlphaFold patterns in {structures_dir} "
                    f"(recursive={search_recursive})."
                )
            out_paths.append(hit)
        return out_paths

    if mode == "colabfold":
        colabfold = shutil.which("colabfold_batch")
        if not colabfold:
            raise RuntimeError(
                "colabfold_batch not found in PATH; install ColabFold or use --structure-mode existing."
            )
        out_dir = structures_dir / f"colabfold_{fasta_path.stem}"
        out_dir.mkdir(parents=True, exist_ok=True)
        # Prepare a clean FASTA with sanitized headers matching seq.name
        prep_fa = out_dir / f"{fasta_path.stem}.sanitized.fasta"
        with open(prep_fa, "w", encoding="utf-8") as g:
            for rec in seqs:
                g.write(f">{rec.name}\n{wrap_seq(rec.seq)}\n")
        # Run colabfold
        cmd = [colabfold, str(prep_fa), str(out_dir)]
        cp = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if cp.returncode != 0:
            raise RuntimeError(
                f"colabfold_batch failed for {fasta_path.name}:\n{cp.stderr}"
            )
        # Collect top-ranked PDB for each sequence (heuristic search)
        for rec in seqs:
            candidates = sorted(out_dir.rglob(f"*{rec.name}*rank*_model*.pdb"))
            if not candidates:
                candidates = sorted(out_dir.rglob(f"{rec.name}.pdb"))
            if not candidates:
                candidates = sorted(out_dir.rglob("*.pdb"))
            if not candidates:
                raise FileNotFoundError(
                    f"No PDB produced by ColabFold for sequence {rec.name}."
                )
            out_paths.append(candidates[0])
        return out_paths

    raise ValueError(f"Unsupported structure mode: {mode}")


# ----------------------------- Gene Ontology uppropagation ----------------
# DF Columns = ["Protein", "interpro", "Score", "GO ID"]
def process_go_uppropagation(
    go_obo: str,
    organism: str,
    aggregate_rows: pd.DataFrame,
):
    """
    Process GO up-propagation for a given set of annotations.

    Args:
        go_obo (str): Path to the Gene Ontology OBO file.
        organism (str): A name for the annotation set (e.g., 'my_organism').
        aggregate_rows (pd.DataFrame): DataFrame with columns ["Protein", "GO ID", "Score"].

    Returns:
        pd.DataFrame: A new DataFrame with up-propagated GO annotations.
    """
    if aggregate_rows.empty:
        return pd.DataFrame(columns=aggregate_rows.columns)

    print(f"[INFO] Aggregate rows shape: {aggregate_rows.shape}")

    go = GeneOntology.GeneOntology(go_obo, verbose=True)
    go.build_structure()

    organism_name = organism if organism else "default_organism"

    print("[INFO] Loading annotations for up-propagation...")
    # Ensure the DataFrame has the required columns for load_annotations
    annotations_to_load = aggregate_rows[["Protein", "GO ID", "Score"]].copy()
    go.load_annotations(annotations_to_load, organism_name)

    print("[INFO] Up-propagating GO terms...")
    go.up_propagate_annotations(organism_name)

    print("[INFO] Retrieving up-propagated annotations...")
    uppropagated_df = go.get_annotations(organism_name)

    # Merge with original data to retain the 'interpro' column
    # We use a left merge to keep all uppropagated terms,
    # and fill missing 'interpro' values (for new terms) with an empty string.
    # We need to handle cases where a (Protein, GO ID) pair might appear multiple times
    # with different 'interpro' values in the original data.
    # To do this, we'll drop duplicates from the mapping table before merging.
    interpro_mapping = aggregate_rows[
        ["Protein", "GO ID", "interpro"]
    ].drop_duplicates()

    merged_df = pd.merge(
        uppropagated_df,
        interpro_mapping,
        on=["Protein", "GO ID"],
        how="left",
    )
    merged_df["interpro"] = merged_df["interpro"].fillna("")

    print(f"[INFO] Up-propagated DataFrame shape: {merged_df.shape}")
    # Reorder columns to match the desired output format
    return merged_df[["Protein", "interpro", "Score", "GO ID"]]


# ----------------------------- Pipeline Orchestration ------------------------


def main():
    # Argument parsing
    ap = argparse.ArgumentParser(
        description="Batch S2F pipeline: FASTA → structure → Foldseek → InterPro → GO."  # noqa
    )
    g_in = ap.add_mutually_exclusive_group(required=True)
    g_in.add_argument(
        "--fastas",
        nargs="+",
        help="Input FASTA files (glob accepted by shell).",  # noqa
    )
    g_in.add_argument(
        "--query-structures",
        nargs="+",
        help="Directly provide query PDB/mmCIF files (skip structure prediction).",  # noqa
    )

    ap.add_argument(
        "--structure-mode",
        choices=["existing", "colabfold"],
        default="existing",
        help="How to obtain structures for FASTA inputs (default: existing). Ignored when --query-structures is used.",  # noqa
    )
    ap.add_argument(
        "--structures-dir",
        type=Path,
        default=Path("structures"),
        help="Directory to look for (or write) structures for FASTA inputs.",
    )

    ap.add_argument(
        "--target-db",
        required=True,
        help="Foldseek target DB base path (from `foldseek createdb`, e.g., /data/targetDB_afsp).",  # noqa
    )
    ap.add_argument(
        "--protein2ipr", required=True, help="protein2ipr.dat (optionally .gz)."  # noqa
    )
    ap.add_argument(
        "--interpro2go",
        default=None,
        help="(Optional) interpro2go (optionally .gz) to include GO terms.",
    )
    ap.add_argument(
        "--alignment-type",
        type=int,
        default=1,
        choices=[0, 1],
        help="Foldseek alignment type: 0=local 3Di+AA, 1=TM-align re-score/rank (default: 1).",  # noqa
    )
    ap.add_argument(
        "--evalue-max",
        type=float,
        default=1e-5,
        help="Max E-value (default: 1e-5).",  # noqa
    )
    ap.add_argument(
        "--min-qcov",
        type=float,
        default=0.70,
        help="Min query coverage (default: 0.70).",
    )
    ap.add_argument(
        "--min-tcov",
        type=float,
        default=0.70,
        help="Min target coverage (default: 0.70).",
    )
    ap.add_argument(
        "--min-avg-tm",
        type=float,
        default=0.50,
        help="Min average TM-score if available (default: 0.50).",
    )
    ap.add_argument(
        "--topk",
        type=int,
        default=50,
        help="Max accepted non-redundant hits to keep per query (default: 50).",  # noqa
    )
    ap.add_argument(
        "--outdir", type=Path, default=Path("results"), help="Output directory."  # noqa
    )

    ap.add_argument(
        "--search-recursive",
        action="store_true",
        help="Search structures_dir recursively for AlphaFold files.",
    )
    ap.add_argument(
        "--prostt5-model",
        default=None,
        help="Enable FASTA→3Di mode (no query structures). "
        "Pass the ProstT5 weights directory (download via `foldseek databases ProstT5 <dir> <tmp>`).",  # noqa
    )
    ap.add_argument(
        "--gpu",
        default=None,
        help="GPU device list for Foldseek (e.g., '0' or '0,1'). "
        "Requires a CUDA-enabled Foldseek build.",
    )
    ap.add_argument(
        "--go-obo",
        default="go.obo",
        help="Gene Ontology OBO file (default: go.obo).",
    )
    ap.add_argument(
        "--gaf",
        default="goa_uniprot_all.gaf.gz",
        help="GOA UniProt GAF file (optionally .gz) for GO annotations (default: goa_uniprot_all.gaf.gz).",  # noqa
    )
    ap.add_argument(
        "--organism",
        default=None,
        help="(Optional) organism name to filter GOA annotations (e.g., 'Homo sapiens').",  # noqa
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    # ---------------------------------------------------------------

    # Load optional InterPro->GO
    print("[INFO] Loading InterPro2GO mappings...")
    ipr2go = parse_interpro2go(args.interpro2go)

    # Prepare list of query structure files
    # FASTA inputs if any
    query_inputs: List[Path] = []  # type: ignore
    all_hits_per_query = defaultdict(list)

    # Structure inputs if any
    query_structures: List[Path] = []  # type: ignore
    query_names: List[str] = []  # type: ignore

    if args.prostt5_model:
        print("[INFO] ProstT5 FASTA→3Di mode enabled.")
        if not args.fastas:
            sys.exit(
                "When --prostt5-model is given, provide FASTA(s) via --fastas."
            )  # noqa
        # Each FASTA file is contains all of the query sequences per organism
        # (to be processed individually later)
        for fa in args.fastas:
            p = Path(fa)
            if not p.exists():
                sys.exit(f"FASTA not found: {p}")
            query_inputs.append(p)
            query_names.append(p.stem)
    else:
        print("[INFO] Preparing query structures...")
        if args.query_structures:
            for p in args.query_structures:
                print(f"[INFO] Using provided structure: {p}")
                pp = Path(p)
                if not pp.exists():
                    sys.exit(f"Query structure not found: {pp}")
                query_structures.append(pp)
                query_names.append(pp.stem)
        else:
            # From FASTAs, ensure structures
            for fasta in args.fastas:
                print(f"[INFO] Processing FASTA: {fasta}")
                fasta_path = Path(fasta)
                if not fasta_path.exists():
                    sys.exit(f"FASTA not found: {fasta_path}")
                try:
                    struct_paths = ensure_structure_for_fasta(
                        fasta_path,
                        args.structure_mode,
                        args.structures_dir,
                        search_recursive=args.search_recursive,
                    )
                except Exception as e:
                    sys.exit(str(e))
                for sp in struct_paths:
                    query_structures.append(sp)
                    query_names.append(sp.stem)

    # 1) Run Foldseek for each query
    # TODO: For now, this process is working only for the ProstT5 mode
    # We may want to extend it to structure queries too.
    # The variable used for the structures is query_structures
    # but it is not used anywhere else.
    start_time = time.time()
    print("[INFO] Running Foldseek")
    with tempfile.TemporaryDirectory(prefix="s2f_foldseek_") as td:
        tmpdir = Path(td)
        all_hits_per_query: Dict[str, List[Hit]] = {}  # type: ignore
        for qpath, qname in zip(query_inputs, query_names):
            try:
                hits = run_foldseek_easy_search(
                    qpath,
                    args.target_db,
                    tmpdir,
                    alignment_type=args.alignment_type,
                    max_hits=max(200, args.topk),
                    prostt5_model=args.prostt5_model,
                    gpu=args.gpu,
                )
            except Exception as e:
                sys.stderr.write(
                    f"[WARN] Foldseek failed for {qpath.name}: {e}\n"
                )  # noqa
                hits = []
            print(f"[INFO] {len(hits)} hits from Foldseek for {qname}.")
            print("[INFO]: Grouping hits by query sequence...")
            # Regroup by per-sequence query ID
            for qid, subhits in group_hits_by_query(hits).items():
                # breakpoint()
                all_hits_per_query[qid] = subhits

    end_time = time.time()
    print(f"[INFO] Foldseek done in {end_time - start_time:.1f} seconds.")
    # 2) Apply selection / thresholds and
    # collect wanted UniProt accessions across all queries
    print("[INFO] Applying selection thresholds...")
    tm_available = args.prostt5_model is None
    accepted_per_query: Dict[str, List[Hit]] = {}  # type: ignore
    wanted_accs: Set[str] = set()  # type: ignore
    for qid, hits in all_hits_per_query.items():
        acc_hits = select_hits(
            hits,
            evalue_max=args.evalue_max,
            min_qcov=args.min_qcov,
            min_tcov=args.min_tcov,
            min_avg_tm=args.min_avg_tm,
            topk=args.topk,
            tm_available=tm_available,
        )
        accepted_per_query[qid] = acc_hits
        for h in acc_hits:
            acc = extract_uniprot_from_id(h.target)
            if acc:
                wanted_accs.add(acc)

    # 3) Stream protein2ipr only once
    start_time = time.time()
    print(
        f"[INFO] Fetching InterPro mappings for {len(wanted_accs)} unique accessions..."  # noqa
    )
    acc2ipr = stream_protein2ipr_for_wanted(args.protein2ipr, wanted_accs)
    end_time = time.time()
    print(f"[INFO] InterPro mappings loaded in {end_time - start_time:.1f} seconds.")
    # 4) Per-query outputs + build aggregate rows
    # (one row per final IPR with support count)
    print("[INFO] Writing foldseek summary...")
    aggregate_rows = []  # type: ignore
    query_id_map = {}
    for qid in accepted_per_query.keys():
        query_id_map[qid] = normalise_query_identifier(qid)

    for qid, acc_hits in accepted_per_query.items():
        protein_id = query_id_map[qid]
        safe = sanitize_name(qid)
        # Evidence rows for TSV
        rows = []
        ipr_support = Counter()
        for h in acc_hits:
            acc = extract_uniprot_from_id(h.target)
            iprs = acc2ipr.get(acc, set()) if acc else set()
            rows.append(
                [
                    h.target,
                    acc or "",
                    f"{h.evalue:.2e}",
                    f"{h.bits:.1f}",
                    f"{h.qcov:.3f}",
                    f"{h.tcov:.3f}",
                    (f"{h.avg_tm:.3f}" if h.avg_tm is not None else ""),
                    str(len(iprs)),
                ]
            )
            for ipr in iprs:
                ipr_support[ipr] += 1

        # Write per-query TSV
        tsv_path = args.outdir / f"{safe}.hits.tsv"
        with open(tsv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                [
                    "target",
                    "uniprot",
                    "evalue",
                    "bits",
                    "qcov",
                    "tcov",
                    "avg_tm",
                    "ipr_count",
                ]
            )
            w.writerows(rows)

        # Write per-query JSON
        json_path = args.outdir / f"{safe}.summary.json"
        summary = {
            "query": safe,
            "protein_id": protein_id,
            "parameters": {
                "evalue_max": args.evalue_max,
                "min_qcov": args.min_qcov,
                "min_tcov": args.min_tcov,
                "min_avg_tm": args.min_avg_tm,
                "topk": args.topk,
                "alignment_type": args.alignment_type,
                "prostt5_model": bool(args.prostt5_model),
            },
            "accepted_hits": [
                {
                    "target": h.target,
                    "uniprot": extract_uniprot_from_id(h.target) or "",
                    "evalue": h.evalue,
                    "bits": h.bits,
                    "qcov": round(h.qcov, 3),
                    "tcov": round(h.tcov, 3),
                    "avg_tm": (
                        round(h.avg_tm, 3) if h.avg_tm is not None else None
                    ),  # noqa
                }
                for h in acc_hits
            ],
            "interpro_support": dict(ipr_support),
        }

        # Include InterPro->GO if available
        if ipr2go:
            summary["interpro_to_go"] = {
                ipr: sorted(ipr2go.get(ipr, set()))
                for ipr in ipr_support.keys()  # noqa
            }
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    start_time = time.time()
    print("[INFO] Building aggregate assignments...")
    for qid, acc_hits in accepted_per_query.items():
        protein_id = query_id_map[qid]
        ipr_support = Counter()
        for h in acc_hits:
            acc = extract_uniprot_from_id(h.target)
            iprs = acc2ipr.get(acc, set()) if acc else set()
            for ipr in iprs:
                ipr_support[ipr] += 1
            # Aggregate table: one row per IPR,
            # and go term with final support count
            for ipr, _ in ipr_support.items():
                if ipr2go:
                    gos = sorted(ipr2go.get(ipr, set()))
                else:
                    gos = []
                for go in gos:
                    aggregate_rows.append([protein_id, ipr, 1, go])

    end_time = time.time()
    print(f"[INFO] Aggregate assignments built in {end_time - start_time:.1f} seconds.")
    # 6) Aggregate assignment table (query, interpro, support_count, go_terms)
    print("[INFO] Writing aggregate assignments...")
    if aggregate_rows:
        print(
            f"[INFO] {len(aggregate_rows)} total assignments (Probably with duplicates)."
        )
        aggregate_rows_df = pd.DataFrame(
            aggregate_rows,
            columns=["Protein", "interpro", "Score", "GO ID"],  # noqa
        )
        # Drop duplicate (Protein, interpro, GO ID) rows
        aggregate_rows_df = aggregate_rows_df.drop_duplicates(
            subset=["Protein", "interpro", "GO ID"]
        )
        print(
            f"[INFO] {len(aggregate_rows_df)} unique assignments after dropping duplicates."  # noqa
        )  # noqa
        aggregate_rows_df["Score"] = 1  # Reset score to 1 for up-propagation
        aggregate_rows_df = process_go_uppropagation(
            args.go_obo, args.organism, aggregate_rows_df
        )
        print(
            f"[INFO] {len(aggregate_rows_df)} assignments after GO up-propagation."  # noqa
        )  # noqa
        # Create the Protein x GO ID matrix from the up-propagated data
        print(
            "[INFO] Creating Protein x GO ID matrix R from up-propagated data."
        )  # noqa
        matrix_df = aggregate_rows_df.copy()
        matrix_df["Score"] = pd.to_numeric(matrix_df["Score"])

        protein_go_matrix = matrix_df.pivot_table(
            index="Protein",
            columns="GO ID",
            values="Score",
            aggfunc="sum",
            fill_value=0,
        )

        print(
            "[INFO] Protein x GO ID matrix R created with shape:",
            protein_go_matrix.shape,
        )
        agg_path = args.outdir / "assignments.tsv"
        # Rename column for output
        aggregate_rows_df.to_csv(
            agg_path,
            sep="\t",
            header=True,
            index=False,
            quoting=csv.QUOTE_NONE,
        )
        # matrix_path = args.outdir / "R.tsv"
        # protein_go_matrix.to_csv(
        #     matrix_path,
        #     sep="\t",
        #     header=True,
        #     index=True,
        #     quoting=csv.QUOTE_NONE,
        # )

    print(f"[OK] Processed {len(query_structures)} queries.")
    print(f"[OK] Outputs in: {args.outdir}")
    print(
        "[INFO] Per-query files: <query>.hits.tsv and <query>.summary.json; aggregate: assignments.tsv"  # noqa
    )


if __name__ == "__main__":

    start_time = time.time()
    main()
    end_time = time.time()
    print(f"[INFO] Total time taken: {end_time - start_time:.2f} seconds")
