import os
from pyfaidx import Fasta

# download genome files from UCSC
chr_names = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY',
]

# create blank directory to store genome files
local_path = "./genome_files/hg38_fasta/"
if os.path.exists(local_path):
    os.system(f"rm -rf {local_path}")
os.system(f"mkdir -p {local_path}")

# download the files
print("Downloading hg38 genome files from UCSC...")
os.system(f"rsync -az --info=progress2 rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ {local_path}")

print("Unzipping and indexing hg38...")
# iterate through files in genome_files, unzip desired, delete undesired, and build index
files = os.listdir(local_path)
for file in files:
    chr_name = file.split(".")[0]
    if chr_name not in chr_names:
        os.remove(os.path.join(local_path, file))
        continue
    # unzip the file
    os.system(f"gunzip {os.path.join(local_path, file)}")
    Fasta(f"{local_path}/{chr_name}.fa", build_index=True)
print("Done!")
