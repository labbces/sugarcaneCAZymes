# Conserved CAZyme Detector & Graph Exporter

This project reads **dbCAN** annotations and **DIAMOND** pairwise similarity results to detect **conserved CAZymes** across species. It writes:

- A **GraphML** network (`<prefix>.conservedCAZymes.graphml`) with proteins as nodes and conserved hits as edges.
- Tab-delimited reports:
  - `<prefix>.conservedCAZymes.txt` (all conserved CAZyme pairs kept)
  - `<prefix>.conservedCAZymes.targetfams.txt` (subset where at least one side matches configured target families/types)

---

## Overview

Pipeline steps:

1. **dbCAN overview** tables (per species) are read to map each protein to its CAZy class/family/subfamily annotations.
2. **DIAMOND** pairwise comparisons are parsed and filtered by identity, e-value, and coverage.
3. Candidate pairs are consolidated; **bidirectionality** is flagged when the reciprocal is found.
4. Pairs are exported as:
   - a **GraphML** network (`graphml`),
   - **TSV** tables (all pairs; target families only).

---

## Requirements

- Python 3.9+
- Biopython (`Bio.SeqIO`)
- Standard library: `argparse`, `collections`, `glob`, `gzip`, `os`, `re`, `xml.etree`, `xml.dom`
- Graph analysis: networkx and igraph and leidenalg

Create a python virtual environment and activate it:

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
```


Install libraries:

```bash
pip install biopython networkx python-igraph leidenalg
```

---

## Directory Layout (example)

```
project/
├─ data/
│  ├─ dbCAN_results/
│  │  ├─ SpeciesA.dbCAN.overview.tsv.gz
│  │  └─ SpeciesB.dbCAN.overview.tsv.gz
│  ├─ diamond.orig/
│  │  ├─ SpeciesIDs.txt
│  │  ├─ SequenceIDs.txt.gz
│  │  ├─ Blast1_2.txt.gz
│  │  └─ Blast2_1.txt.gz
│  └─ proteins/
│     ├─ SpeciesA.faa.gz
│     ├─ SpeciesA.idx        # created by the script (SeqIO.index_db)
│     └─ SpeciesB.faa.gz
├─ script.py
└─ README.md
```

---

## Input Files

### 1) dbCAN Overview Tables
**Path:** `data/dbCAN_results/*.dbCAN.overview.tsv.gz`  
**Format:** gzip, tab-separated; rows where the **penultimate** column equals `'3'` are kept (high-confidence hits).  
**Species name:** taken from the filename prefix (substring before the first `.`).  
**Output mapping:**  
`{ species_name: { gene_id: "CAZy_annotation|..." } }`

> Example row (kept): `gene123\t...\t3\tGH5|GH5_21`

---

### 2) DIAMOND Results
**Path:** `data/diamond.orig/BlastX_Y.txt.gz`  
**Filename regex:** `Blast(\d+)_(\d+).txt.gz` (X and Y are **species IDs**)  
**Format:** gzip, **outfmt 6 default** columns:
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```
Self vs. self (`X == Y`) are skipped.

---

### 3) Species IDs
**Path:** `data/diamond.orig/SpeciesIDs.txt`  
**Format:** one mapping per line:
```
<int_species_id>:<species_name_or_basename>
```
The `.faa` suffix is removed if present. The resulting species names **must** match the dbCAN species keys and the protein FASTA basenames.

> Example:
> ```
> 1:SpeciesA.faa
> 2:SpeciesB
> ```

---

### 4) Sequence IDs
**Path:** `data/diamond.orig/SequenceIDs.txt.gz`  
**Format:** gzip, one mapping per line:
```
<spid_serial>:<seqname possibly with description>
```
- `spid_serial` looks like `1_000123`.
- Only the **first token** of `<seqname>` (up to first space) is used.
- Kept only if the sequence name exists in dbCAN mapping for that species.

---

### 5) Protein FASTA files
**Path:** `data/proteins/<species>.faa.gz`  
**Format:** gzip-compressed FASTA (one file per species).  
The script creates on-disk indexes (`<species>.idx`) via `Bio.SeqIO.index_db`.

> **Note:** Random access on gz can be slow for very large files; consider BGZF or uncompressed FASTA if performance is an issue.

---

## Output Files

- **GraphML:** `<prefix>.conservedCAZymes.graphml`  
  - Node ID: protein sequence name (FASTA ID)
  - Node attrs: `species`, `seq_len`, `dbcan.subfamilies`, `dbcan.families`, `dbcan.classes`, `dbcan.target_type`
  - Edge attrs: `pident`, `evalue`, `qalgnlen`, `qalgnper`, `salgnlen`, `salgnper`, `bidirectional`, `node1.sp.name`, `node2.sp.name`
  - Default `edgedefault="undirected"` (reciprocals are deduplicated)

- **All conserved pairs (TSV):** `<prefix>.conservedCAZymes.txt`  
  Header:
  ```
  #node1.sp.name  node1.seq.name  node1.seq.len  node1.dbcan.subfamilies  node1.dbcan.families  node1.dbcan.classes  node1.dbcan.target_type  node2.sp.name  node2.seq.name  node2.seq.len  node2.dbcan.subfamilies  node2.dbcan.families  node2.dbcan.classes  node2.dbcan.target_type  pident  evalue  qalgnlen  qalgnper  salgnlen  salgnper  bidirectional
  ```

- **Target families subset (TSV):** `<prefix>.conservedCAZymes.targetfams.txt`  
  Same columns as above; includes rows where **either** side has `dbcan.target_type != 'N/A'`.

---

## Filtering Criteria

A DIAMOND hit is kept if **all** are true:

- `pident >= 80`
- `evalue < 1e-5`
- `query` and `subject` coverage ≥ **50%**  
  (computed from `(qend - qstart + 1)/qlen` and `(send - sstart + 1)/slen`)
- **Bidirectional** is marked `True` if the reverse pair was seen earlier.

---

## Target Family Configuration

Curated in-script via `targetCAZyFamilies`:

- Groups (e.g., **Endo-xilanases**, **Arabinofuranosidases**, …)
- Each group lists CAZy families and optional **subfamilies** to include.
  - `subfamilies = None` means **entire family** qualifies.
- Transformed into `family_to_targets` (family → {group → config}) for fast lookup.

> Example snippet:
> ```python
> targetCAZyFamilies = {
>   'Endo-xilanases': {
>       "GH5":  {"class": "GH", "subfamilies": {21, 34, 35}},
>       "GH10": {"class": "GH", "subfamilies": None},  # include entire family
>       ...
>   },
>   ...
> }
> ```
> The script computes a per-protein `dbcan.target_type` by checking whether any of its subfamilies (or whole family) matches a configured group.

---

## Usage

```bash
python script.py -p OUT_PREFIX [-V 0|1|2|3]
```

- `-p / --prefix` **(required):** Output prefix, e.g., `results/run1`
- `-V / --verbose` *(optional)*: 0 (silent), 1 (normal), 2 (debug), 3 (trace). Default: 1

**Examples:**

```bash
# Normal run
python script.py -p results/cazy_run1

# Debug mode
python script.py -p results/cazy_run1 -V 2
```

---

## Graph Schema (quick reference)

**Node attributes**
- `species` *(str)*  
- `seq_len` *(int)*
- `dbcan.subfamilies` *(str or None; comma-joined if multiple)*
- `dbcan.families` *(str or None; comma-joined if multiple)*
- `dbcan.classes` *(str or None; comma-joined if multiple)*
- `dbcan.target_type` *(str; comma-joined group names or `'N/A'`)*

**Edge attributes**
- `pident` *(double)*
- `evalue` *(double)*
- `qalgnlen` *(int)*, `qalgnper` *(double %)*
- `salgnlen` *(int)*, `salgnper` *(double %)*
- `bidirectional` *(boolean)*
- `node1.sp.name`, `node2.sp.name` *(string; convenience)*

---

## Performance Notes

- `SeqIO.index_db` creates disk indexes (`.idx`) for faster repeated access.
- Very large gz FASTAs can be I/O-bound; prefer **BGZF** or uncompressed FASTA if indexing becomes slow.
- The script skips DIAMOND pairs that aren’t present in the dbCAN-filtered set to reduce lookups.

---

## Troubleshooting

- **Species in `SpeciesIDs.txt` not found in dbCAN results**  
  Ensure the species basename (after removing `.faa`) matches the prefix of the dbCAN file.  
  Example: `1:SpeciesA.faa` → must have `data/dbCAN_results/SpeciesA.dbCAN.overview.tsv.gz`.

- **Missing protein FASTA**  
  Check `data/proteins/<species>.faa.gz` exists for every species used in `SpeciesIDs.txt`.

- **No rows in outputs**  
  Loosen thresholds (identity/coverage/e-value) or verify your dbCAN filters (penultimate column must be `'3'`).

- **Indexing WARNING already present**  
  Safe to ignore; reuses the existing `.idx` files.

---

## License

