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
    # 3 Trace  TRACE    Fine-grained steps. per-TE messages, etc
    if verbose >= level:
        print(msg)

def read_dbcan_tables(folder_path):
    """
    Read dbCAN overview files (pattern: SPP.dbCAN.overview.tsv.gz) and return a mapping of species to selected hits.

    For each file:
      - species_name is taken from the filename prefix (text before the first dot).
      - lines are parsed as tab-separated fields; empty lines are ignored.
      - only keep rows with at least two columns and with the penultimate column equal to '3' (high-confidence hits).
      - for each kept row map the first column (gene/protein ID) to the last column (CAZyme family/annotation).

    Returns:
        dict: { species_name (str): { gene_id (str): annotation (str), ... }, ... }
        Example: { "SPP": { "geneX": "GH5", "geneY": "CE1" }, "Other": { ... } }

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
    return dbcan_dict

def to_graph_structs(keep_dbcan_results, undirected=True):
    """
    Build:
      nodes: id -> attrs  (node id = sequence name)
      edges: [(u, v, attrs), ...]
    Dedupes reciprocal edges if undirected=True.
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
    if isinstance(value, bool):  return "boolean"
    if isinstance(value, int):   return "int"
    if isinstance(value, float): return "double"
    return "string"

def _collect_attr_schema(nodes_dict, edges_list):
    node_keys, edge_keys = {}, {}
    for _, attrs in nodes_dict.items():
        for k, v in (attrs or {}).items():
            node_keys.setdefault(k, _gml_type(v))
    for _, _, attrs in edges_list:
        for k, v in (attrs or {}).items():
            edge_keys.setdefault(k, _gml_type(v))
    return node_keys, edge_keys

def write_graphml(nodes, edges, path, directed=False, default_edge_weight=None):
    NS = {"g": "http://graphml.graphdrawing.org/xmlns"}
    ET.register_namespace("", NS["g"])
    graphml = ET.Element(ET.QName(NS["g"], "graphml"))
    node_schema, edge_schema = _collect_attr_schema(nodes, edges)
    if default_edge_weight is not None and "weight" not in edge_schema:
        edge_schema["weight"] = _gml_type(default_edge_weight)
    key_id = 0
    key_index = {}
    def add_key(name, domain, gtype, default=None):
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
                if family_to_targets[ufam][target_group]['subfamilies'] is None or any(subfam in subfamilies_list for subfam in map(str, family_to_targets[ufam][target_group]['subfamilies'])):
                    ufam_type.add(target_group)
    ufam_type_str = 'N/A'
    if ufam_type:
        ufam_type_str = ','.join(sorted(ufam_type))
    # print(ufam_type_str)
    return classes, sub_families, families, ufam_type_str

def select_diamond_results(diamond_folder='data/diamond', species_ids_file='SpeciesIDs.txt', sequences_ids_file='SequenceIDs.txt.gz', sequences_folder='data/proteins/', dbcan_res=None, family_to_targets=None, verbose=1):
    """
    Read gzip tabular DIAMOND results named BlastX_Y.txt.gz from diamond_folder.
    Skip files where X == Y. Map integer IDs X/Y to species names using species_ids_file.

    select CAZymes conserved across species based on DIAMOND results and dbCAN annotations.
    """
    
    log(f"Reading DIAMOND results from: {diamond_folder} using species IDs from: {species_ids_file} and sequence IDs from: {sequences_ids_file}", 1, verbose)

    keep_dbcan_results = defaultdict(dict)

    # read species id -> name mapping
    # just keep species that are in dbcan_res
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

    # read sequence id -> name mapping
    # just keep sequences that are in dbcan_res
    id_to_seqname = defaultdict(dict)  # especie -> {seqid: seqname}
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
    # read DIAMOND results
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
        
        # Build disk-based index for protein sequences for both species (x_name and y_name)
        # Assumes sequence files are named as "{species}.faa.gz" in a folder, e.g., "data/proteins/"

        seq_indexes = {}

        for sp in [x_name, y_name]:
            seq_file = os.path.join(sequences_folder, f"{sp}.faa.gz") #sequence files must be compressed using bgzip for random access
            if os.path.exists(seq_file):
            # Create an index file for each species using Bio.SeqIO.index_db
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

        idx1=SeqIO.index_db(seq_indexes[x_name]) 
        idx2=SeqIO.index_db(seq_indexes[y_name]) 
        log(f'Reading DIAMOND results for: {x_name} vs {y_name}',1, verbose)

        with gzip.open(file_path, 'rt') as fh:
            #Default fields in diamond output:
            #qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                if len(cols) < 12:
                    continue
                sp1,seqser1=cols[0].split("_")
                sp2,seqser2=cols[1].split("_")
                # print(f'{file_path}\t{line}')
                # print(f'Species 1: {sp1} {id_to_spname[int(sp1)]}, Sequence 1: {cols[0]}')
                # print(f'Species 2: {sp2} {id_to_spname[int(sp2)]}, Sequence 2: {cols[1]}')
                if cols[0] in id_to_seqname[id_to_spname[int(sp1)]] and cols[1] in id_to_seqname[id_to_spname[int(sp2)]]:
                    seqname1=id_to_seqname[id_to_spname[int(sp1)]][cols[0]]
                    seqname2=id_to_seqname[id_to_spname[int(sp2)]][cols[1]]
                    # print(f'seqname1: {seqname1}, seqname2: {seqname2}')
                    seq_record1 = idx1[seqname1]
                    seq_record2 = idx2[seqname2]
                    seq_len1 = len(seq_record1.seq)
                    seq_len2 = len(seq_record2.seq)
                    qalgnlen=int(cols[7]) - int(cols[6]) + 1
                    salgnlen=int(cols[9]) - int(cols[8]) + 1
                    qalgnper=qalgnlen/seq_len1 * 100
                    salgnper=salgnlen/seq_len2 * 100
                    spseqname1 = f'{id_to_spname[int(sp1)]}__{seqname1}'
                    spseqname2 = f'{id_to_spname[int(sp2)]}__{seqname2}'
                    #filter based on alignment metrics
                    #keeping hits with at least 80% identity, e-value < 1e-5, and at least 50% alignment coverage in both sequences
                    bidirectional = False
                    if float(cols[2]) >= 80 and float(cols[10]) < 1e-5 and qalgnper >= 50 and salgnper >= 50:
                        if (spseqname2 in keep_dbcan_results) and (spseqname1 in keep_dbcan_results[spseqname2]):
                                bidirectional = True
                                keep_dbcan_results[spseqname2][spseqname1]['bidirectional'] = bidirectional
                                continue
                        subfamilies1_list=dbcan_res[id_to_spname[int(sp1)]][seqname1].split('|')
                        classes1, sub_families1, families1, ufam1_type_str =produce_families_strs(subfamilies1_list,family_to_targets)
                        subfamilies2_list=dbcan_res[id_to_spname[int(sp2)]][seqname2].split('|')
                        classes2, sub_families2, families2, ufam2_type_str =produce_families_strs(subfamilies2_list,family_to_targets)

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
    return keep_dbcan_results

args = parse_arguments()
targetCAZyFamilies = {
    #notes for Igor: GH43 subfamilies 7 and 16 are both in Endo-Xylanases and Arabinofuranosidases
    #There are no subfamilies for GH38 in CAZY
    #GH 43 only have 40 subfamilie sin CAZY, there is no subfamily 177
    'Endo-xilanases' :{
        "GH5":  {"class": "GH", "subfamilies": {21, 34, 35}},
        "GH10": {"class": "GH", "subfamilies": None},  # include entire family
        "GH11": {"class": "GH", "subfamilies": None},
        "GH8": {"class": "GH", "subfamilies": None},
        "GH38": {"class": "GH", "subfamilies": {7, 8}}, 
        "GH43": {"class": "GH", "subfamilies": {7, 16, 177}},
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
    'alpha-Galactosidases':{
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
        family_to_targets[fam][group] = entry


dbcan_res=read_dbcan_tables(folder_path='data/dbCAN_results/')
# print("Read dbCAN results:", dbcan_res.keys())
dbcan_res = select_diamond_results(diamond_folder='data/diamond.orig', species_ids_file='SpeciesIDs.txt', sequences_ids_file='SequenceIDs.txt.gz', sequences_folder='data/proteins/', dbcan_res=dbcan_res, family_to_targets=family_to_targets,verbose=args.verbose)
#Conserved CAZymes out file
nodes, edges =to_graph_structs(dbcan_res, undirected=True)
write_graphml(nodes, edges, f"{args.prefix}.conservedCAZymes.graphml", directed=False)
conserved_cazymes_file = f"{args.prefix}.conservedCAZymes.txt"
conserved_cazymes_targetfams_file = f"{args.prefix}.conservedCAZymes.targetfams.txt"
out_conserved_cazymes = open(conserved_cazymes_file, 'wt')
out_conserved_cazymes_targetfams = open(conserved_cazymes_targetfams_file, 'wt')

log(f"Writing conserved CAZymes to: {conserved_cazymes_file}", 1, args.verbose)
print('#node1.sp.name\tnode1.seq.name\tnode1.seq.len\tnode1.dbcan.subfamilies\tnode1.dbcan.families\tnode1.dbcan.classes\tnode1.dbcan.target_type\tnode2.sp.name\tnode2.seq.name\tnode2.seq.len\tnode2.dbcan.subfamilies\tnode2.dbcan.families\tnode2.dbcan.classes\tnode2.dbcan.target_type\tpident\tevalue\tqalgnlen\tqalgnper\tsalgnlen\tsalgnper\tbidirectional', file=out_conserved_cazymes)
for node1 in dbcan_res:
    for node2 in dbcan_res[node1]:
            res_str=(f'{dbcan_res[node1][node2]["node1.sp.name"]}\t'
                  f'{dbcan_res[node1][node2]["node1.seq.name"]}\t'
                  f'{dbcan_res[node1][node2]["node1.seq.len"]}\t'
                  f'{dbcan_res[node1][node2]["node1.dbcan.subfamilies"]}\t'
                  f'{dbcan_res[node1][node2]["node1.dbcan.families"]}\t'
                  f'{dbcan_res[node1][node2]["node1.dbcan.classes"]}\t'
                  f'{dbcan_res[node1][node2]["node1.dbcan.target_type"]}\t'
                  f'{dbcan_res[node1][node2]["node2.sp.name"]}\t'
                  f'{dbcan_res[node1][node2]["node2.seq.name"]}\t'
                  f'{dbcan_res[node1][node2]["node2.seq.len"]}\t'
                  f'{dbcan_res[node1][node2]["node2.dbcan.subfamilies"]}\t'
                  f'{dbcan_res[node1][node2]["node2.dbcan.families"]}\t'
                  f'{dbcan_res[node1][node2]["node2.dbcan.classes"]}\t'
                  f'{dbcan_res[node1][node2]["node2.dbcan.target_type"]}\t'
                  f'{dbcan_res[node1][node2]["pident"]}\t'
                  f'{dbcan_res[node1][node2]["evalue"]}\t'
                  f'{dbcan_res[node1][node2]["qalgnlen"]}\t'
                  f'{dbcan_res[node1][node2]["qalgnper"]:.1f}\t'
                  f'{dbcan_res[node1][node2]["salgnlen"]}\t'
                  f'{dbcan_res[node1][node2]["salgnper"]:.1f}\t'
                  f'{dbcan_res[node1][node2]["bidirectional"]}')
            print(res_str, file=out_conserved_cazymes)
            if dbcan_res[node1][node2]["node1.dbcan.target_type"] != 'N/A' or dbcan_res[node1][node2]["node2.dbcan.target_type"] != 'N/A':
                print(res_str, file=out_conserved_cazymes_targetfams)