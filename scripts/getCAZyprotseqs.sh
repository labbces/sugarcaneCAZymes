# 0) Install Entrez Direct
# https://www.ncbi.nlm.nih.gov/books/NBK179288/
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=$HOME/edirect:$PATH

# (Recommended) Set an NCBI API key to increase rate limits (free)
export NCBI_API_KEY="d6f7ec46daa0dce15d77dc1c4c4dc61ec109"

# 1) Split your 3M list into chunks of 1000 (tune as needed)
mkdir -p chunks out
split -l 1000 cazy_data_20250707.acc.nr chunks/acc_

# 2) Fetch in parallel (keep concurrency sane: e.g., 4–8 jobs)
# Each chunk posts 1000 IDs, then retrieves FASTA.
# Retries included to handle transient failures.
parallel -j 6 --joblog fetch.log '
  for i in {1..5}; do
    if epost -db protein -input {} \
      | efetch -format fasta > out/{/.}.fa 2>> out/{/.}.err; then
      break
    else
      sleep $((RANDOM%5+2))
    fi
  done
' ::: chunks/acc_*

# 3) Concatenate and compress
cat out/*.fa > cazy_proteins.fasta
bgzip -@ 8 cazy_proteins.fasta

# 4) Quick sanity check
grep -c "^>" cazy_proteins.fasta.gz
