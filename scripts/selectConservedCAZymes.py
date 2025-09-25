# -*- coding: utf-8 -*-
"""
Pipeline to detect conserved CAZymes across species using:
  1) dbCAN functional annotation tables (per-species)
  2) DIAMOND all-vs-all tabular alignments (per species pair)
  3) FASTA protein sequences per species

Inputs (by function):
- read_dbcan_tables(folder_path):
    Expects files: <species>.dbCAN.overview.tsv.gz
    Format (tsv, gz): keep rows where the penultimate column == '3' (high-confidence).
    Maps: gene/protein ID -> final dbCAN family annotation (last column; pipe-separated subfamilies).

- select_diamond_results():
    diamond_folder:
        - SpeciesIDs.txt
            Lines like: 1:SpeciesA.faa
            (integer id):(species fasta name, .faa suffix will be stripped)
        - SequenceIDs.txt.gz
            Lines like: 1_000001:SeqName1 some extra description
            (speciesID_seqSerial):(sequence identifier and description)
            Only the first token after ':' is taken as the sequence name.
        - BlastX_Y.txt.gz
            DIAMOND tabular (outfmt 6) per species pair X vs Y.
            Default 12 cols: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    sequences_folder:
        - One gzipped FASTA per species: <species>.faa.gz
        - An on-disk index will be created: <species>.idx (Bio.SeqIO.index_db)

Outputs:
- GraphML file: <prefix>.conservedCAZymes.graphml
- Tabular files:
    * <prefix>.conservedCAZymes.txt              (all conserved)
    * <prefix>.conservedCAZymes.targetfams.txt   (subset where the CAZy families match your target groups)
"""
from collections import defaultdict
import os
import glob
import gzip
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from Bio import SeqIO
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generates a distilled TE file from genome-specific runs of TE annotation in several species')
    parser.add_argument('--verbose', '-V', type=int, default=1, help='Set verbosity level (0 = silent, 1 = normal, 2 = debug). Default is 1.')
    parser.add_argument('--prefix', '-p', required=True, type=str, help='Prefix used to create output files.')
    return parser.parse_args()

def log (msg, level=1,  verbose=0):
    #Default verbose level is 1 (normal)
    #Verbose levels should be one of:
    # 0 Silent CRITICAL Only critical errors
    # 1 Normal INFO     Key steps, progress messages
    # 2 Debug  DEBUG    Extra info: file paths, filtered seqs
    # 3 Trace  TRACE    Fine-grained steps
    if verbose >= level:
        print(msg)

def read_dbcan_tables(folder_path,summary_file=None, verbose=1):
    """
    Read dbCAN "overview" files and select high-confidence hits.

    Input:
      folder_path/ <species>.dbCAN.overview.tsv.gz
        - gzip compressed TSV
        - Keep rows with penultimate column == '3' (high confidence)
        - Map: first column (gene/protein ID) -> last column (CAZy annotated families/subfamilies)
          Note: last column may be pipe-separated subfamilies (e.g., GH5_21|GH5_34|...).

    Returns:
      dict: { species_name: { gene_id: annotation_str, ... }, ... }
      Example: { "SpeciesA": { "geneX": "GH5_21|GH5_34", "geneY": "CE1" }, ... }
    """
    dbcan_dict = {}
    pattern = os.path.join(folder_path, "*.dbCAN.overview.tsv.gz")
    for file_path in glob.glob(pattern):
        species_name = os.path.basename(file_path).split('.')[0]
        with gzip.open(file_path, 'rt') as f:
            rows = [line.strip().split('\t') for line in f if line.strip()]
            filtered_rows = [row for row in rows if len(row) >= 2 and row[-2] == '3']
            rows = {row[0]: row[-1] for row in filtered_rows}
        dbcan_dict[species_name] = rows
    with open(summary_file, 'wt') as out_summary:
        print("Species\t#Number Proteins with dbCAN annotation", file=out_summary)
        for sp in sorted(dbcan_dict.keys()):
            print(f"{sp}\t{len(dbcan_dict[sp])}", file=out_summary)
    log(f"Read dbCAN annotations for {len(dbcan_dict)} species from {folder_path}", 1, verbose)
    return dbcan_dict

def to_graph_structs(keep_dbcan_results, undirected=True):
    """
    Convert conserved hits into graph data structures.

    Input:
      keep_dbcan_results: nested dict keyed by "Species__SeqName"
        {
          "SpA__Seq1": {
            "SpB__Seq2": {
               'node1.sp.name': 'SpA', 'node1.seq.name': 'Seq1', ...
               'node2.sp.name': 'SpB', 'node2.seq.name': 'Seq2', ...
               'pident': float, 'evalue': float, ...
            },
            ...
          },
          ...
        }

    Returns:
      nodes: dict { node_id(str) -> attrs(dict) }
             node_id is the plain sequence name (not "Species__Seq")
      edges: list [ (u, v, attrs_dict), ... ]
             If undirected=True, reciprocal edges are deduplicated by (sorted(u,v)).
    """
    nodes = {}
    edges = []
    seen = set()  # for deduping edges when undirected

    for seq1, inner in keep_dbcan_results.items():
        for seq2, rec in inner.items():
            # Node IDs: use the sequence names (unique and handy)
            n1 = rec['node1.seq.name']
            n2 = rec['node2.seq.name']

            # Node attributes (species, length, dbCAN summaries)
            if n1 not in nodes:
                nodes[n1] = {
                    'species': rec['node1.sp.name'],
                    'seq_len': rec['node1.seq.len'],
                    'dbcan.subfamilies': rec['node1.dbcan.subfamilies'],
                    'dbcan.families': rec['node1.dbcan.families'],
                    'dbcan.classes': rec['node1.dbcan.classes'],
                    'dbcan.target_type': rec['node1.dbcan.target_type'],
                }
            if n2 not in nodes:
                nodes[n2] = {
                    'species': rec['node2.sp.name'],
                    'seq_len': rec['node2.seq.len'],
                    'dbcan.subfamilies': rec['node2.dbcan.subfamilies'],
                    'dbcan.families': rec['node2.dbcan.families'],
                    'dbcan.classes': rec['node2.dbcan.classes'],
                    'dbcan.target_type': rec['node2.dbcan.target_type'],
                }

            # Edge attributes (carry alignment metrics + bidirectional flag)
            eattrs = {
                'pident': rec['pident'],
                'evalue': rec['evalue'],
                'qalgnlen': rec['qalgnlen'],
                'qalgnper': rec['qalgnper'],
                'salgnlen': rec['salgnlen'],
                'salgnper': rec['salgnper'],
                'bidirectional': rec['bidirectional'],
                # Optional helpful context:
                'node1.sp.name': rec['node1.sp.name'],
                'node2.sp.name': rec['node2.sp.name'],
            }

            if undirected:
                key = tuple(sorted((n1, n2)))
                if key in seen:
                    # If you want, merge/upgrade attributes here (e.g., set bidirectional=True)
                    continue
                seen.add(key)
                edges.append((key[0], key[1], eattrs))
            else:
                edges.append((n1, n2, eattrs))

    return nodes, edges

def _gml_type(value):
    """Infer GraphML attribute types."""
    if isinstance(value, bool):  return "boolean"
    if isinstance(value, int):   return "int"
    if isinstance(value, float): return "double"
    return "string"

def _collect_attr_schema(nodes_dict, edges_list):
    """
    Inspect node/edge attributes to define GraphML key schema
    (attribute names + types).
    """
    node_keys, edge_keys = {}, {}
    for _, attrs in nodes_dict.items():
        for k, v in (attrs or {}).items():
            node_keys.setdefault(k, _gml_type(v))
    for _, _, attrs in edges_list:
        for k, v in (attrs or {}).items():
            edge_keys.setdefault(k, _gml_type(v))
    return node_keys, edge_keys

def write_graphml(nodes, edges, path, directed=False, default_edge_weight=None):
    """
    Serialize nodes/edges into GraphML.

    Args:
      nodes: dict {node_id: {attr: value}}
      edges: list of (u, v, attrs_dict)
      path: output .graphml filename (str)
      directed: bool (GraphML 'edgedefault')
      default_edge_weight: if provided, inject 'weight' for edges missing it
    """
    NS = {"g": "http://graphml.graphdrawing.org/xmlns"}
    ET.register_namespace("", NS["g"])
    graphml = ET.Element(ET.QName(NS["g"], "graphml"))
    node_schema, edge_schema = _collect_attr_schema(nodes, edges)
    if default_edge_weight is not None and "weight" not in edge_schema:
        edge_schema["weight"] = _gml_type(default_edge_weight)
    key_id = 0
    key_index = {}
    def add_key(name, domain, gtype, default=None):
        """
        Register <key> in GraphML header and store mapping (domain, name) -> key id.
        """
        nonlocal key_id
        kid = f"k{key_id}"; key_id += 1
        key_el = ET.SubElement(
            graphml, ET.QName(NS["g"], "key"),
            id=kid, **{"for": domain}, attr_name=name, attr_type=gtype
        )
        # xml.etree does not allow hyphens in attrib names unless using dict unpack trick above
        key_el.set("attr.name", name)
        key_el.set("attr.type", gtype)
        if default is not None:
            d = ET.SubElement(key_el, ET.QName(NS["g"], "default"))
            d.text = str(default)
        key_index[(domain, name)] = kid
    for name, gtype in sorted(node_schema.items()):
        add_key(name, "node", gtype)
    for name, gtype in sorted(edge_schema.items()):
        add_key(name, "edge", gtype)
    graph = ET.SubElement(graphml, ET.QName(NS["g"], "graph"),
                          edgedefault="directed" if directed else "undirected")
    for nid, attrs in nodes.items():
        n = ET.SubElement(graph, ET.QName(NS["g"], "node"), id=str(nid))
        for k, v in (attrs or {}).items():
            kid = key_index.get(("node", k))
            if kid is None:
                add_key(k, "node", _gml_type(v))
                kid = key_index[("node", k)]
            d = ET.SubElement(n, ET.QName(NS["g"], "data"), key=kid)
            d.text = str(v)
    for i, (u, v, *rest) in enumerate(edges):
        attrs = rest[0] if rest else {}
        e = ET.SubElement(graph, ET.QName(NS["g"], "edge"),
                          id=f"e{i}", source=str(u), target=str(v))
        if default_edge_weight is not None and "weight" not in (attrs or {}):
            kid = key_index.get(("edge", "weight"))
            d = ET.SubElement(e, ET.QName(NS["g"], "data"), key=kid); d.text = str(default_edge_weight)
        for k, v in (attrs or {}).items():
            kid = key_index.get(("edge", k))
            if kid is None:
                add_key(k, "edge", _gml_type(v))
                kid = key_index[("edge", k)]
            d = ET.SubElement(e, ET.QName(NS["g"], "data"), key=kid)
            d.text = str(v)
    rough = ET.tostring(graphml, encoding="utf-8")
    pretty = minidom.parseString(rough).toprettyxml(indent="  ", encoding="utf-8")
    with open(path, "wb") as fh:
        fh.write(pretty)

def produce_families_strs(subfamilies_list=None, family_to_targets=None):
    """
    Summarize dbCAN subfamilies into:
      - classes   (e.g., 'GH', 'CE'), possibly comma-joined if mixed
      - sub_families (string with comma-joined GHx_y terms, or single, or None)
      - families  (e.g., 'GH5', 'CE1'), possibly comma-joined if mixed
      - ufam_type_str (comma-joined target group labels from 'family_to_targets')

    Args:
      subfamilies_list: list like ['GH5_21','GH5_34'] (parsed from dbCAN string split by '|')
      family_to_targets: dict { 'GH5': { 'Endo-xilanases': {'class': 'GH', 'subfamilies': {21,34,35}}, ...}, ... }

    Returns:
      (classes, sub_families, families, ufam_type_str)
    """
    unique_classes = set()
    unique_families = set()
    for fam in subfamilies_list:
        match = re.match(r'^(([A-Za-z]+)\d+)_?', fam)
        if match:
            unique_classes.add(match.group(2))
            unique_families.add(match.group(1))
    if len(unique_classes) > 1:
        classes = ','.join(sorted(unique_classes))
    else:
        classes = next(iter(unique_classes)) if unique_classes else None
    if len(subfamilies_list) > 1:
        sub_families = ','.join(sorted(subfamilies_list))
    else:
        sub_families = next(iter(subfamilies_list)) if subfamilies_list else None
    if len(unique_families) > 1:
        families = ','.join(sorted(unique_families))
    else:
        families = next(iter(unique_families)) if unique_families else None
    
    ufam_type= set()
    for ufam in unique_families:
        if ufam in family_to_targets:
            for target_group in family_to_targets[ufam]:
                if family_to_targets[ufam][target_group]['subfamilies'] is None:
                    ufam_type.add(target_group)
                else:
                    for subfam in subfamilies_list:
                        if subfam.startswith(ufam+'_'):
                            subfam_num = subfam.split('_')[1]
                            if subfam_num in map(str, family_to_targets[ufam][target_group]['subfamilies']):
                                ufam_type.add(target_group)

    ufam_type_str = 'N/A'
    if ufam_type:
        ufam_type_str = ','.join(sorted(ufam_type))
    # print(ufam_type_str)
    return classes, sub_families, families, ufam_type_str

def select_diamond_results(diamond_folder='data/diamond', 
                           species_ids_file='SpeciesIDs.txt', 
                           sequences_ids_file='SequenceIDs.txt.gz', 
                           sequences_folder='data/proteins/', 
                           dbcan_res=None, 
                           family_to_targets=None, 
                           verbose=1):
    """
    Integrate DIAMOND results + dbCAN to find conserved CAZymes.

    Input files (located under 'diamond_folder'):
      - SpeciesIDs.txt (text)
          Lines: "<int_id>:<species_name>.faa"
          Example: "1:Zea_mays.faa"
          Note: ".faa" is stripped to get species_name.
          Must match species present in dbCAN overview files already read.

      - SequenceIDs.txt.gz (gzip)
          Lines: "<spid>_<serial>:<seqname> [optional description]"
          Example: "1_000001:Zm00001d034567.1 protein..."
          Only the first token after ':' (split by space) is used as the sequence name.

      - BlastX_Y.txt.gz (gzip) files:
          DIAMOND outfmt 6 with standard columns:
          qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
          Where qseqid/sseqid are "<spid>_<serial>"

    Protein sequences (sequences_folder):
      - One gzipped FASTA per species: "<species_name>.faa.gz"
      - Index file "<species_name>.idx" is created/used by Bio.SeqIO.index_db for random access.

    Filtering criteria for conserved hits:
      - pident >= 80
      - evalue < 1e-5
      - query/subject coverage >= 50% of their respective protein lengths
      - bidirectional flag set if reciprocal hit already present in the reverse file.

    Returns:
      keep_dbcan_results (nested dict) as consumed by 'to_graph_structs'.
    """
    
    log(f"Reading DIAMOND results from: {diamond_folder} using species IDs from: {species_ids_file} and sequence IDs from: {sequences_ids_file}", 1, verbose)

    keep_dbcan_results = defaultdict(dict)
    count_kepth_hits=defaultdict(lambda: defaultdict(int))
    sp803280_conserved_hits = defaultdict(
            lambda: defaultdict(
                lambda: {
                    "queryseqinfo": {},
                    "hits": defaultdict(lambda: defaultdict(dict))  # hits[sp2][seqname2] -> dict
                }))
    species_dict={}

    # ------------------------------------------------------------
    # 1) Read SpeciesIDs.txt: map int ID -> species name
    #    Only species present in dbCAN results are kept.
    # ------------------------------------------------------------
    id_to_spname = {}
    try:
        species_ids_path = os.path.join(diamond_folder, species_ids_file)
        with open(species_ids_path, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split(":")
                if len(parts) < 2:
                    continue
                try:
                    sid = int(parts[0])
                except ValueError:
                    continue
                parts[1] = parts[1].replace('.faa', '')
                # Check if species name exists in dbcan_dict
                species_name = parts[1].strip()
                if species_name not in dbcan_res.keys():
                    raise ValueError(f"Species name '{species_name}' from SpeciesIDs.txt not found in dbcan_dict. Stop processing.")
                id_to_spname[sid] = species_name
    except FileNotFoundError:
        print(f"Error: Species IDs file '{species_ids_file}' not found. Stopping processing.")
        raise

    # ------------------------------------------------------------
    # 2) Read SequenceIDs.txt.gz: map "<spid>_<serial>" -> "<seqname>"
    #    Only sequences that appear in dbCAN for that species are kept.
    # ------------------------------------------------------------

    id_to_seqname = defaultdict(dict) # species -> { '<spid>_<serial>': 'SeqName' }
    try:
        sequences_ids_path = os.path.join(diamond_folder, sequences_ids_file)
        with gzip.open(sequences_ids_path, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                seqid, seqname = line.split(":",1)
                spid,seqserial=seqid.split("_",1)
                spname=id_to_spname[int(spid)]
                seqid = seqid.strip()
                seqname = seqname.strip().split(' ')[0]  # take only the first part of the sequence name
                if spname is None:
                    print(f"Warning: Species ID {spid} not found in id_to_spname. Skipping sequence ID {seqid}.")
                    exit(1)
                elif spname in dbcan_res.keys():
                    if seqname in dbcan_res[spname].keys():
                        # print(f'Sequence ID: {seqid}, Sequence Name: {seqname}, Species ID: {spid}, Species Name: {id_to_spname[int(spid)]}')
                        id_to_seqname[spname][seqid]=seqname
    except FileNotFoundError:
        print(f"Error: Sequence IDs file '{sequences_ids_file}' not found. Stopping processing.")
        raise

    for i in id_to_spname:
        log(f'C:{i}\t{id_to_spname[int(i)]}\t{len(id_to_seqname[id_to_spname[int(i)]].keys())}',2, verbose=1)

    # ------------------------------------------------------------
    # 3) Iterate DIAMOND BlastX_Y.txt.gz files and collect conserved hits
    # ------------------------------------------------------------
    pattern = os.path.join(diamond_folder, "Blast*_*.txt.gz")
    for file_path in glob.glob(pattern):
        fname = os.path.basename(file_path)
        m = re.match(r'Blast(\d+)_(\d+)\.txt\.gz$', fname)
        if not m:
            continue
        x_id, y_id = int(m.group(1)), int(m.group(2))
        if x_id == y_id:
            continue
        x_name = id_to_spname.get(x_id, str(x_id))
        y_name = id_to_spname.get(y_id, str(y_id))
        species_dict[x_name]=1
        species_dict[y_name]=1
        
        # ------------------------------------------------------------
        # Build (or reuse) on-disk FASTA indexes for both species
        # Expect "<species>.faa.gz" files in sequences_folder
        # ------------------------------------------------------------

        seq_indexes = {}

        for sp in [x_name, y_name]:
            seq_file = os.path.join(sequences_folder, f"{sp}.faa.gz") #sequence files must be compressed using bgzip for random access
            if os.path.exists(seq_file):
                # Create a persistent index database for random access to sequences
                index_db_path = os.path.join(sequences_folder, f"{sp}.idx")
                log(f'{index_db_path}',2,verbose)
                if not os.path.exists(index_db_path):
                    SeqIO.index_db(index_db_path, seq_file, "fasta")
                    seq_indexes[sp] = index_db_path
                else:
                    print(f"Warning: Sequence index file for species '{sp}' already present at {index_db_path}.")
                    seq_indexes[sp] = index_db_path
            else:
                print(f"Warning: Sequence file for species '{sp}' not found at {seq_file}. Stopping now. Make sure all protein files are present in the {sequences_folder} folder.")
                exit(1)
 
        # Open the two species indexes
        idx1=SeqIO.index_db(seq_indexes[x_name]) 
        idx2=SeqIO.index_db(seq_indexes[y_name]) 

        log(f'Reading DIAMOND results for: {x_name} vs {y_name}',1, verbose)

        # Read DIAMOND lines and apply filters
        with gzip.open(file_path, 'rt') as fh:
            # outfmt 6 (12 default columns):
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                if len(cols) < 12:
                    continue
                sp1,seqser1=cols[0].split("_")
                sp2,seqser2=cols[1].split("_")

                # Keep only if both qseqid and sseqid are known in dbCAN mapping
                if cols[0] in id_to_seqname[id_to_spname[int(sp1)]]and cols[1] in id_to_seqname[id_to_spname[int(sp2)]]:
                    seqname1=id_to_seqname[id_to_spname[int(sp1)]][cols[0]]
                    seqname2=id_to_seqname[id_to_spname[int(sp2)]][cols[1]]

                    # Protein lengths via indexed FASTA (random access)
                    seq_record1 = idx1[seqname1]
                    seq_record2 = idx2[seqname2]
                    seq_len1 = len(seq_record1.seq)
                    seq_len2 = len(seq_record2.seq)

                    # Alignment spans and coverage (%)
                    qalgnlen=int(cols[7]) - int(cols[6]) + 1
                    salgnlen=int(cols[9]) - int(cols[8]) + 1
                    qalgnper=qalgnlen/seq_len1 * 100
                    salgnper=salgnlen/seq_len2 * 100

                    # Compose compound keys: "Species__SeqName"
                    spseqname1 = f'{id_to_spname[int(sp1)]}__{seqname1}'
                    spseqname2 = f'{id_to_spname[int(sp2)]}__{seqname2}'

                    # filter based on alignment metrics
                    # keeping hits with at least 80% identity, e-value < 1e-5, and at least 50% alignment coverage in both sequences
                    bidirectional = False
                    if float(cols[2]) >= 80 and float(cols[10]) < 1e-5 and qalgnper >= 50 and salgnper >= 50:
                        # Parse dbCAN subfamilies (pipe-separated string) for both nodes
                        subfamilies1_list=dbcan_res[id_to_spname[int(sp1)]][seqname1].split('|')
                        classes1, sub_families1, families1, ufam1_type_str =produce_families_strs(subfamilies1_list,family_to_targets)
                        subfamilies2_list=dbcan_res[id_to_spname[int(sp2)]][seqname2].split('|')
                        classes2, sub_families2, families2, ufam2_type_str =produce_families_strs(subfamilies2_list,family_to_targets)

                        if id_to_spname[int(sp1)] == 'SP803280':
                            entry = sp803280_conserved_hits[id_to_spname[int(sp1)]][seqname1]
                            entry["queryseqinfo"].update({
                                'sp1.dbcan.subfamilies': sub_families1,
                                'sp1.dbcan.families': families1,
                                'sp1.dbcan.classes': classes1,
                                'sp1.dbcan.target_type': ufam1_type_str,
                            })
                            # print(entry)
                            entry["hits"][id_to_spname[int(sp2)] ][seqname2] = {
                                'pident': float(cols[2]),
                                'evalue': float(cols[10]),
                                'qalgnlen': qalgnlen,
                                'salgnlen': salgnlen,
                                'qalgnper': qalgnper,
                                'salgnper': salgnper,
                                'sp2.dbcan.subfamilies': sub_families2,
                                'sp2.dbcan.families':    families2,
                                'sp2.dbcan.classes':     classes2,
                                'sp2.dbcan.target_type': ufam2_type_str,
                            }
                        # If we already stored the reverse edge, mark it as bidirectional and skip adding the duplicate
                        if (spseqname2 in keep_dbcan_results) and (spseqname1 in keep_dbcan_results[spseqname2]):
                                bidirectional = True
                                keep_dbcan_results[spseqname2][spseqname1]['bidirectional'] = bidirectional
                                continue
                        
                        count_kepth_hits[id_to_spname[int(sp1)]][id_to_spname[int(sp2)]]+=1
                        # Store record
                        keep_dbcan_results[spseqname1][spseqname2] = {
                            'node1.sp.name': id_to_spname[int(sp1)],
                            'node1.seq.name': seqname1,
                            'node1.seq.len': seq_len1,
                            'node1.dbcan.subfamilies': sub_families1,
                            'node1.dbcan.families': families1,
                            'node1.dbcan.classes': classes1,
                            'node1.dbcan.target_type': ufam1_type_str,
                            'node2.sp.name': id_to_spname[int(sp2)],
                            'node2.seq.name': seqname2,
                            'node2.seq.len': seq_len2,
                            'node2.dbcan.subfamilies': sub_families2,
                            'node2.dbcan.families': families2,
                            'node2.dbcan.classes': classes2,
                            'node2.dbcan.target_type': ufam2_type_str,
                            'pident': float(cols[2]),
                            'evalue': float(cols[10]),
                            'qalgnlen': qalgnlen,
                            'qalgnper': qalgnper,
                            'salgnlen': salgnlen,
                            'salgnper': salgnper,
                            'bidirectional': bidirectional
                        }
                        # print(f'{id_to_spname[int(sp1)]}\t{seqname1}\t{dbcan_res[id_to_spname[int(sp1)]][seqname1]}\t{id_to_spname[int(sp1)]}\t{seq_len1}\t{seqname2}\t{dbcan_res[id_to_spname[int(sp2)]][seqname2]}\t{cols[2]}\t{cols[10]}\t{seq_len2}\t{qalgnlen}\t{qalgnper:.1f}\t{salgnlen}\t{salgnper:.1f}', file=out_conserved_cazymes)
    return keep_dbcan_results, count_kepth_hits, sp803280_conserved_hits, species_dict

# ======================
# MAIN EXECUTION
# ======================
args = parse_arguments()
targetCAZyFamilies = {
    'Endo-xilanases' :{
        "GH5":  {"class": "GH", "subfamilies": {21, 34, 35}},
        "GH10": {"class": "GH", "subfamilies": None},  # include entire family
        "GH11": {"class": "GH", "subfamilies": None},
        "GH8": {"class": "GH", "subfamilies": None},
        "GH30": {"class": "GH", "subfamilies": {7, 8}}, 
        "GH43": {"class": "GH", "subfamilies": {7, 16, 17}}, # also in Arabinofuranosidases
        "GH98": {"class": "GH", "subfamilies": None},
        "GH141": {"class": "GH", "subfamilies": None},
        "CE1":  {"class": "CE", "subfamilies": None}, 
    },
    'Arabinofuranosidases':{
        "GH43": {"class": "GH", "subfamilies": {2, 7, 9, 10, 12, 16, 17, 18, 19, 20, 21, 22, 23, 26, 27, 29, 33, 34, 35, 36}},
        "GH51": {"class": "GH", "subfamilies": {1, 2}},
        "GH54": {"class": "GH", "subfamilies": { None }},
        "GH159": {"class": "GH", "subfamilies": { None }},
    },
    'Glucuronidases':{
        "GH2": {"class": "GH", "subfamilies": { None }},
        "GH79": {"class": "GH", "subfamilies": { None }},
    },
    'Acetil-xilano estereases':{
        "CE1": {"class": "CE", "subfamilies": { None }},
        "CE2": {"class": "CE", "subfamilies": { None }},
        "CE3": {"class": "CE", "subfamilies": { None }},
        "CE6": {"class": "CE", "subfamilies": { None }},
        "CE7": {"class": "CE", "subfamilies": { None }},
    },
    'Endo-mananases':{
        "GH5": {"class": "GH", "subfamilies": {4, 7, 8, 10, 17, 18, 19, 25, 31, 36, 40, 41, 55}},
        "GH26": {"class": "GH", "subfamilies": { None }},
        "GH113": {"class": "GH", "subfamilies": { None }},
        "GH134": {"class": "GH", "subfamilies": { None }},
    },
    'Alpha-Galactosidases':{
        "GH27": {"class": "GH", "subfamilies": { None }},
        "GH36": {"class": "GH", "subfamilies": { None }},
        "GH57": {"class": "GH", "subfamilies": { None }},
        "GH97": {"class": "GH", "subfamilies": { None }},
        "GH4": {"class": "GH", "subfamilies": { None }},
        "GH110": {"class": "GH", "subfamilies": { None }},
    },
}
# Transform targetCAZyFamilies to use family (e.g., "GH5") as main key
family_to_targets = {}
for group, families in targetCAZyFamilies.items():
    for fam, info in families.items():
        if fam not in family_to_targets:
            family_to_targets[fam] = {}
        entry = info.copy()
        subs = entry.get('subfamilies')
        # If it's {None} or empty set → treat as None (whole family)
        if subs == {None} or subs == set():
            entry['subfamilies'] = None
            
        family_to_targets[fam][group] = entry

# ----------------------------------------
# Load dbCAN overview results for all species
# ----------------------------------------
summary_file=f'{args.prefix}.dbCAN.summary.txt'
dbcan_res=read_dbcan_tables(folder_path='data/dbCAN_results/',summary_file=summary_file, verbose=args.verbose)

# ----------------------------------------
# Read DIAMOND results, filter, and consolidate conserved CAZymes
# ----------------------------------------
conserv_dbcan_res, count_kepth_hits, sp803280_kepth_hits, species = select_diamond_results(diamond_folder='data/diamond.orig', 
                                   species_ids_file='SpeciesIDs.txt', 
                                   sequences_ids_file='SequenceIDs.txt.gz', 
                                   sequences_folder='data/proteins/', 
                                   dbcan_res=dbcan_res, 
                                   family_to_targets=family_to_targets,
                                   verbose=args.verbose)

# ----------------------------------------
# Build graph and write GraphML
# ----------------------------------------
nodes, edges =to_graph_structs(conserv_dbcan_res, undirected=True)
write_graphml(nodes, edges, f"{args.prefix}.conservedCAZymes.graphml", directed=False)

# ----------------------------------------
# Write tabular outputs
# ----------------------------------------
conserved_cazymes_file = f"{args.prefix}.conservedCAZymes.txt"
conserved_cazymes_targetfams_file = f"{args.prefix}.conservedCAZymes.targetfams.txt"

log(f"Writing summary to: {summary_file}", 1, args.verbose)
with open(summary_file, 'at') as out_summary:
    out_summary.write('\n##########################\n#Summary of conserved CAZymes between species pairs\n')
    out_summary.write('\nSpecies1\tSpecies2\t#Number Conserved CAZymes\n')
    for sp1 in count_kepth_hits:
        for sp2 in count_kepth_hits[sp1]:
            out_summary.write(f'{sp1}\t{sp2}\t{count_kepth_hits[sp1][sp2]}\n')

log(f"Writing conserved CAZymes to: {conserved_cazymes_file}", 1, args.verbose)

#write results table for SP803280 conserved CAZymes
sp803280_conserved_cazymes_file = f"{args.prefix}.SP803280.conservedCAZymes.txt"
log(f"Writing SP803280 conserved CAZymes to: {sp803280_conserved_cazymes_file}", 1, args.verbose)
sorted_species = sorted(species)

with open(sp803280_conserved_cazymes_file, "wt") as out_sp803280_conserved_cazymes:
    # Header: sp1, seqid_sp1, then one column per species
    header = "sp1\tseqid_sp1(fam)\ttarget\t" + "\t".join(sorted_species) + "\n"
    out_sp803280_conserved_cazymes.write(header)

    for sp1, seqdict in sp803280_kepth_hits.items():
        for seqid_sp1, entry in seqdict.items():
            q_subs  = entry["queryseqinfo"].get("sp1.dbcan.subfamilies", "")
            q_ttype = entry["queryseqinfo"].get("sp1.dbcan.target_type", "")
            seq1_annot = f"{seqid_sp1}({q_subs})"
            row_cells = [sp1, seq1_annot, q_ttype]

            for sp2 in sorted_species:
                if sp2 in entry["hits"]:
                    annotated = []
                    # deterministic order of target seq IDs
                    for seqid_sp2, hitinfo in sorted(entry["hits"][sp2].items()):
                        subs = hitinfo.get("sp2.dbcan.subfamilies", "")
                        # fams = hitinfo.get("sp2.dbcan.families", "")
                        # clss = hitinfo.get("sp2.dbcan.classes", "")
                        # annotated.append(f"{seqid_sp2}({subs};{fams};{clss})")
                        annotated.append(f"{seqid_sp2}({subs})")
                    cell = ",".join(annotated)
                else:
                    cell = ""  # keep column alignment

                row_cells.append(cell)

            out_sp803280_conserved_cazymes.write("\t".join(row_cells) + "\n")



# Emit records for all conserved pairs; also write a filtered file for target families
with open(conserved_cazymes_file, 'wt') as out_conserved_cazymes, open(conserved_cazymes_targetfams_file, 'wt') as out_conserved_cazymes_targetfams:
    # Header line for main table
    out_conserved_cazymes.write('#node1.sp.name\tnode1.seq.name\tnode1.seq.len\tnode1.dbcan.subfamilies\tnode1.dbcan.families\tnode1.dbcan.classes\tnode1.dbcan.target_type\tnode2.sp.name\tnode2.seq.name\tnode2.seq.len\tnode2.dbcan.subfamilies\tnode2.dbcan.families\tnode2.dbcan.classes\tnode2.dbcan.target_type\tpident\tevalue\tqalgnlen\tqalgnper\tsalgnlen\tsalgnper\tbidirectional\n')
    out_conserved_cazymes_targetfams.write('#node1.sp.name\tnode1.seq.name\tnode1.seq.len\tnode1.dbcan.subfamilies\tnode1.dbcan.families\tnode1.dbcan.classes\tnode1.dbcan.target_type\tnode2.sp.name\tnode2.seq.name\tnode2.seq.len\tnode2.dbcan.subfamilies\tnode2.dbcan.families\tnode2.dbcan.classes\tnode2.dbcan.target_type\tpident\tevalue\tqalgnlen\tqalgnper\tsalgnlen\tsalgnper\tbidirectional\n')
    for node1 in conserv_dbcan_res:
        for node2 in conserv_dbcan_res[node1]:
                res_str=(f'{conserv_dbcan_res[node1][node2]["node1.sp.name"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node1.seq.name"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node1.seq.len"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node1.dbcan.subfamilies"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node1.dbcan.families"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node1.dbcan.classes"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node1.dbcan.target_type"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.sp.name"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.seq.name"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.seq.len"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.dbcan.subfamilies"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.dbcan.families"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.dbcan.classes"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["node2.dbcan.target_type"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["pident"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["evalue"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["qalgnlen"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["qalgnper"]:.1f}\t'
                    f'{conserv_dbcan_res[node1][node2]["salgnlen"]}\t'
                    f'{conserv_dbcan_res[node1][node2]["salgnper"]:.1f}\t'
                    f'{conserv_dbcan_res[node1][node2]["bidirectional"]}')
                out_conserved_cazymes.write(res_str)
                # If either side matches a configured target group, also write to the target-families file
                if conserv_dbcan_res[node1][node2]["node1.dbcan.target_type"] != 'N/A' or conserv_dbcan_res[node1][node2]["node2.dbcan.target_type"] != 'N/A':
                    out_conserved_cazymes_targetfams.write(res_str)