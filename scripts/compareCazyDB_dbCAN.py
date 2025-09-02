import glob
import gzip

dbcan_results = {}

# Use the provided path pattern
for filename in glob.glob("/data/dmrp/cazy/cazy_proteins_CAZymes_*/overview.tsv"):
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            # Check if last but one column equals 3
            if parts[-2] == "3":
                key = parts[0]
                value = parts[-1]
                dbcan_results[key] = value



cazy_classification = {}

with gzip.open("cazy_data_20250707.txt.gz", "rt") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) != 5:
            continue
        if parts[-1] == "ncbi":
            key = parts[-2]
            value = parts[0]
            cazy_classification[key] = value

for key in dbcan_results:
    if key in cazy_classification:
        print(f"{key}\t{dbcan_results[key]}\t{cazy_classification[key]}")