# Instructions for `foldseek_interpro2go.py`

## Summary

The `foldseek_interpro2go.py` script is an automated pipeline for protein function annotation. It takes a protein sequence (FASTA) or structure (PDB/CIF) as input, finds structurally similar proteins using `foldseek`, and transfers their functional annotations—specifically InterPro domains and Gene Ontology (GO) terms—to the original query protein.

This "Structure-to-Function" (S2F) approach is powerful for annotating proteins that lack clear sequence homology to known proteins but may share a similar fold.

## Key Features

-   **Flexible Input**: Accepts protein sequences (FASTA) or 3D structures (PDB/CIF).
-   **Structure Prediction**: Can automatically predict structures from sequences using ColabFold if a structure is not provided.
-   **Structural Search**: Leverages the speed and sensitivity of `foldseek` to find structural homologs.
-   **Functional Mapping**: Maps structural hits to InterPro domains and subsequently to GO terms.
-   **Rich Output**: Generates detailed per-query reports (TSV and JSON) and a final summary table of all functional assignments.
-   **Sequence-Only Mode**: Supports a `ProstT5` mode to run `foldseek` on sequences without needing full 3D structures.

## Dependencies

### Software
1.  **`foldseek`**: Must be installed and available in your system's `PATH`.
2.  **`colabfold_batch`** (Optional): Required only if you use the `--structure-mode colabfold` option.

### Data Files
1.  **Foldseek Target Database**: A database created with `foldseek createdb` (e.g., against the AlphaFold DB).
2.  **`protein2ipr.dat`**: The mapping file from UniProt accessions to InterPro domains. You can download it from the InterPro FTP.
3.  **`interpro2go`** (Optional): The mapping file from InterPro domains to GO terms. Required if you want GO term annotations.
4.  **ProstT5 Weights** (Optional): Required for the sequence-only (`--prostt5-model`) mode. Download them using `foldseek`:
    ```bash
    foldseek databases ProstT5 ./prostt5_weights ./tmp
    ```

## Usage

The script is run from the command line. Below are the main arguments and some examples.

### Core Arguments

| Argument               | Description                                                                                             |
| :--------------------- | :------------------------------------------------------------------------------------------------------ |
| `--fastas`             | One or more input FASTA files.                                                                          |
| `--query-structures`   | One or more input structure files (PDB/CIF). Use this if you already have the structures.               |
| `--target-db`          | **Required**. Path to the Foldseek target database.                                                     |
| `--protein2ipr`        | **Required**. Path to the `protein2ipr.dat` file (can be gzipped).                                      |
| `--interpro2go`        | (Optional) Path to the `interpro2go` file to enable GO term mapping.                                    |
| `--outdir`             | Directory to save the output files (default: `results`).                                                |
| `--structure-mode`     | How to get structures for FASTA inputs: `existing` (default) or `colabfold`.                            |
| `--structures-dir`     | Directory to find existing structures or save predicted ones (default: `structures`).                   |
| `--evalue-max`         | Maximum E-value for a Foldseek hit to be considered significant (default: `1e-5`).                      |
| `--topk`               | Maximum number of non-redundant hits to process per query (default: `50`).                              |

### Example Commands

**1. Using a pre-existing structure file:**

```bash
python foldseek_interpro2go.py \
    --query-structures /path/to/my_protein.pdb \
    --target-db /path/to/foldseek_db/afdb \
    --protein2ipr /path/to/data/protein2ipr.dat.gz \
    --interpro2go /path/to/data/interpro2go \
    --outdir ./my_protein_results
```

**2. Using a FASTA file and predicting the structure with ColabFold:**

```bash
python foldseek_interpro2go.py \
    --fastas /path/to/my_sequence.fasta \
    --structure-mode colabfold \
    --structures-dir ./predicted_structures \
    --target-db /path/to/foldseek_db/afdb \
    --protein2ipr /path/to/data/protein2ipr.dat.gz \
    --outdir ./my_sequence_results
```

**3. Using a FASTA file and searching for an existing AlphaFold structure:**

```bash
python foldseek_interpro2go.py \
    --fastas /path/to/my_sequence.fasta \
    --structure-mode existing \
    --structures-dir /path/to/alphafold_database/ \
    --search-recursive \
    --target-db /path/to/foldseek_db/afdb \
    --protein2ipr /path/to/data/protein2ipr.dat.gz \
    --outdir ./my_sequence_results
```

**4. Using a FASTA file in sequence-only mode (ProstT5):**

```bash
python foldseek_interpro2go.py \
    --fastas /path/to/my_sequence.fasta \
    --prostt5-model /path/to/prostt5_weights/ \
    --target-db /path/to/foldseek_db/afdb \
    --protein2ipr /path/to/data/protein2ipr.dat.gz \
    --outdir ./my_sequence_results_prostt5
```

## Output Files

The script will generate the following files in the specified output directory:

1.  **`<query_name>.hits.tsv`**: A tab-separated file listing each significant structural hit, its similarity scores (E-value, coverage, TM-score), and the number of InterPro domains associated with it.

2.  **`<query_name>.summary.json`**: A JSON file containing:
    - The parameters used for the run.
    - A list of all accepted hits and their properties.
    - `interpro_support`: A dictionary mapping each found InterPro domain to its "support count" (how many hits contained that domain).
    - `interpro_to_go`: A dictionary mapping the found InterPro domains to their corresponding GO terms (if `--interpro2go` was used).

3.  **`assignments.tsv`**: An aggregate table summarizing the results for all queries. Each row contains a query, an assigned InterPro domain, its support count, and the corresponding GO terms. This file is useful for a high-level overview of the functional annotations.
