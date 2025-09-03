import gzip

dbcan_results = {}

with gzip.open("cazy_proteins_dbcan.tsv.gz", "rt") as f:
    for line in f:
        if line.startswith("Gene ID"):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        # Check if last but one column equals 3
        if parts[-2] == "3":
            SeqID = parts[0]
            family = parts[-1]
            dbcan_results[SeqID] = family

cazy_classification = {}

with gzip.open("cazy_data_20250707.txt.gz", "rt") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) != 5:
            continue
        if parts[-1] == "ncbi":
            seqID = parts[-2]
            family = parts[0]            
            if seqID not in cazy_classification:
                cazy_classification[seqID] = {}
            cazy_classification[seqID][family] = True


combo_dict = {}

with open("dbcan_vs_cazy_families_perseq.tsv", "a") as out_f:
    for key in dbcan_results:
        if key in cazy_classification:
            comb_key = (dbcan_results[key], cazy_classification[key])
            combo_dict.setdefault(dbcan_results[key], set()).update(cazy_classification[key].keys())
            out_f.write(f"{key}\t{dbcan_results[key]}\t{','.join(cazy_classification[key].keys())}\n")
        # print(f"{key}\t{dbcan_results[key]}\t{', '.join(cazy_classification[key].keys())}")


combo_dict = dict(sorted(combo_dict.items()))
with open("dbcan_family_vs_cazy_families.tsv", "w") as summary_f:
    for dbcan_family, cazy_families in combo_dict.items():
        cazy_families_str = ",".join(sorted(cazy_families))
        summary_f.write(f"{dbcan_family}\t{cazy_families_str}\n")
        # print(f"{dbcan_family}\t{cazy_families_str}")