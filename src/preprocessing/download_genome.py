import os
from pyfaidx import Fasta

import subprocess

def run_rsync(os_command, max_attempts=5):
    # run up to five times if rsync command fails
    for i in range(max_attempts):
        result = subprocess.run(os_command, shell=True, check=False)
        if result.returncode == 0:
            break
        if i == max_attempts - 1:
            raise RuntimeError(f"Rsync failed to pull files from UCSC after {max_attempts} attempts. Alert Developer")
        print(f"Attempt {i + 1} failed. Retrying...")
        
    return

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


# download the hg38 chr files
print("Downloading hg38 genome files from UCSC...")
rsync_command = f"rsync -az --info=progress2 --include='chr?.fa.gz' --include='chr??.fa.gz' --exclude='*' rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ {local_path}" 
run_rsync(rsync_command)

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


files = os.listdir(local_path)
fa_files = [f for f in files if f.endswith(".fa")]
fai_files = [f for f in files if f.endswith(".fai")]
if len(fa_files) != 24 or len(fai_files) != 24:
    print(f"Warning: hg38 files not properly downloaded. {len(fa_files)} fa files and {len(fai_files)} fai files found.")
    raise RuntimeError("hg38 files not properly downloaded. Alert Developer")


print("Done!")
