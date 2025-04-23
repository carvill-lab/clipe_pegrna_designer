import os
import subprocess

def run_rsync(os_command, max_attempts=5):
    # run up to five times if rsync command fails
    for i in range(max_attempts):
        result = subprocess.run(os_command, shell=True, check=False)
        if result.returncode == 0:
            break
        if i == max_attempts - 1:
            raise RuntimeError(f"Rsync failed to pull files from NCBI after {max_attempts} attempts. Alert Developer")
        print(f"Attempt {i + 1} failed. Retrying...")
        
    return


print("Downloading latest ClinVar file")
# create blank directory to store clinvar files
local_path = "./genome_files/clinvar/"
if os.path.exists(local_path):
    os.system(f"rm -rf {local_path}")
os.system(f"mkdir -p {local_path}")

# pull the latest clinvar file
rsync_command = f"rsync -az --info=progress2 --exclude='*papu*' --include='clinvar_*.vcf.gz' --include='clinvar_*.tbi' --exclude='*' rsync://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ {local_path}"
run_rsync(rsync_command)
# pull file name
files = os.listdir(local_path)
vcf_files = [f for f in files if f.endswith(".vcf.gz")]
tbi_files = [f for f in files if f.endswith(".tbi")]
if len(vcf_files) != 1 or len(tbi_files) != 1:
    print(f"Warning: clinvar files not properly downloaded. {len(vcf_files)} vcf files and {len(tbi_files)} tbi files found.")
    raise RuntimeError("Clinvar files not properly downloaded. Alert Developer")


# --------

print("Done!")
