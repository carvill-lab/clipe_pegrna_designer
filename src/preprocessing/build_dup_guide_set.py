import concurrent.futures
import gzip
import os
import pickle
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

# TODO limit to 100 bp before and after coding regions
def preprocess_pams(path_str):
    genome_fastas = list(Path(path_str).rglob("*.fa"))
    unique_spacers = set()
    bad_spacers = set()

    # iterate through all chromosomes
    genome_progress = tqdm(genome_fastas)
    for fasta_path in genome_progress:
        genome_progress.set_description(fasta_path.stem)
        fasta_str = str(SeqIO.read(fasta_path, "fasta").seq).upper()
        # iterate through the fasta string looking for potential PAMs
        chr_progress = tqdm(range(20, len(fasta_str)-23), leave=False, desc="searching for pams", unit="b", unit_scale=True, unit_divisor=1000000)
        for i in chr_progress:
            potential_pam = fasta_str[i:i+3]
            if potential_pam[1:] == "GG":
                spacer = fasta_str[i-20:i]
            elif potential_pam[:-1] == "CC":
                spacer = str(Seq(fasta_str[i+3:i+23]).reverse_complement())
            else:
                continue
            
            # spacer found -- process it
            if "N" in spacer or spacer in bad_spacers:
                continue

            if spacer in unique_spacers:
                unique_spacers.remove(spacer)
                bad_spacers.add(spacer)
            else:
                unique_spacers.add(spacer)

    print(f"num unique spacers: {len(unique_spacers)}")
    print(f"non-unique, bad spacers: {len(bad_spacers)}")
    return bad_spacers

    
def main():
    print("iterating through chromosome fasta files to find unique spacers")
    print("this may take a while, and will use >20GB of memory")
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        non_uniques = executor.submit(preprocess_pams, './genome_files/hg38_fasta/').result()
    
    print("exporting non-unique spacers with pickle")
    with open('./genome_files/bad_guides.pkl',  'wb') as handle:
        pickle.dump(non_uniques, handle)
    print("zipping pkl file")
    os.system("gzip -9 ./genome_files/bad_guides.pkl")

if __name__ == "__main__":
    main()

# def load_spacers():
#     with gzip.open('./genome_files/bad_guides.pkl.gz',  'rb') as handle:
#         bad_spacers = pickle.load(handle)
#     print(f"non-unique, bad spacers: {len(bad_spacers)}")

#     # print the first 10 bad spacers
#     print("first 10 bad spacers:")
#     for i, spacer in enumerate(bad_spacers):
#         if i == 10:
#             break
#         print(spacer)

#load_spacers()