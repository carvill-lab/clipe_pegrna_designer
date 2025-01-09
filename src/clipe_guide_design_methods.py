import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from pathlib import Path
pd.options.mode.copy_on_write = True

class clipe_expt:
    def __init__ (self, transcript_name, clinvar_path, gnomad_path, pbs_len, rtt_len, num_windows, disrupt_pam, design_strategy, inclusion_types, allele_count_min=5, prog_bar=None):
        # input values
        self.pbs_len = int(pbs_len)
        self.rtt_len = int(rtt_len)
        self.num_windows = int(num_windows)
        self.disrupt_pam = bool(disrupt_pam)
        self.inclusion_types = list(inclusion_types)
        self.allele_count_min = allele_count_min
        self.prog_bar = prog_bar

        # process input files
        if self.prog_bar:
            self.prog_bar.set(1, detail="Processing inputs")
        self.var_df, self.coding_strand = self.process_input_files(clinvar_path, gnomad_path, transcript_name)
        self.set_variant_locations()

        # pull in correct fasta
        fasta_dir = str(Path(__file__).parent) +  "/genome_files/"
        self.set_ref_fasta(fasta_dir)

        # build dataframes with path, benign, and vus
        self.build_desired_dfs()

        # set desired design strategy
        if design_strategy == "pos_control":
            self.desired_vars = pd.concat([self.pathogenic_vars, self.benign_vars])
        elif design_strategy == "vus":
            self.desired_vars = self.VUS_vars
        else:
            print("Design strategy not recognized")
            raise ValueError("Design strategy not recognized")
        
    def run_guide_design(self):
        # find the densest windows
        if self.prog_bar:
            self.prog_bar.set(3, detail="Finding variant dense regions")
        top_windows = self.find_variant_dense_windows(self.desired_vars)
        include_ptcs = False
        # based on inclusion types, add additional vars to windows
        vars_to_include = self.desired_vars
        if "GNOMAD" in self.inclusion_types and not self.gnomad_vars.empty:
            vars_to_include = pd.concat([vars_to_include, self.gnomad_vars])
        if "BLB" in self.inclusion_types:
            vars_to_include = pd.concat([vars_to_include, self.benign_vars])
        if "PLP" in self.inclusion_types:
            vars_to_include = pd.concat([vars_to_include, self.pathogenic_vars])
        if "PTC" in self.inclusion_types:
            include_ptcs = True
        # remove duplicates
        vars_to_include = vars_to_include.drop_duplicates(subset="var_id")

        if self.prog_bar:
            self.prog_bar.set(4, detail="Designing pegRNAs")
        # design pegRNAs for each window
        final_peg_df = pd.DataFrame()
        guide_screening_df = pd.DataFrame()
        for x, window in enumerate(top_windows):
            window_df = vars_to_include[(vars_to_include["pos"] >= window['rtt_start']) & (vars_to_include["pos"] <= window['rtt_end'])]
            peg_df = self.design_pegrnas(window['rtt_start'], window['rtt_end'], window['peg_strand'], window_df, stop_codons=include_ptcs) 
            peg_df["editing_window"] = x + 1
            final_peg_df = pd.concat([final_peg_df, peg_df])

            temp_screen_df = self.design_pegrnas(window['rtt_start'], window['rtt_end'], window['peg_strand'], window_df, guide_screening_mode=True)
            temp_screen_df["editing_window"] = x + 1
            temp_screen_df['var_id'] = f'window_{x+1}_' + temp_screen_df['var_id']
            guide_screening_df = pd.concat([guide_screening_df, temp_screen_df])

        idt_df = self.clipe_prepare_oligo_df(final_peg_df)
        final_peg_df = pd.merge(final_peg_df, idt_df, on="var_id", how="left")
        
        dfs_to_return = []
        for pegrna_df in [final_peg_df, guide_screening_df]:
            pegrna_df = pegrna_df.drop(columns=["ref_seq", "alt_seq", "read_frame_pos", "index", ])
            pegrna_df = pegrna_df[["editing_window"] + [col for col in pegrna_df.columns if col != "editing_window"]]         
            dfs_to_return.append(pegrna_df)

        return dfs_to_return[0], dfs_to_return[1], top_windows
    

    def set_ref_fasta(self, fasta_dir):
        # Get the unique chromosomes from the variant dataframe
        unique_chromosomes = self.var_df['chr'].unique()
        if len(unique_chromosomes) > 1:
            print("Multiple chromosomes detected in the input file.")
            raise ValueError("Multiple chromosomes detected in the input file.")
        chrom = f'chr{unique_chromosomes[0]}'
        if self.prog_bar:
            self.prog_bar.set(2, detail=f"Loading hg38 {chrom}")
        # Filter the REF_FASTA dictionary to include only the relevant chromosomes
        path = fasta_dir + f"{chrom}.fa.gz"
        self.ref_fasta = str(SeqIO.read(gzip.open(path, "rt"), "fasta").seq).upper()


    def set_variant_locations(self):
        regions = self.var_df.agg({"pos": ["min", "max"]})
        self.vars_start = regions.loc["min", "pos"]
        self.vars_end = regions.loc["max", "pos"]


    def process_input_files(self, clinvar_csv_path, gnomad_csv_path, transcript_name):
        # input validation
        if not clinvar_csv_path.endswith(".txt") and not clinvar_csv_path.endswith(".tsv"):
            raise ValueError("Clinvar file must be a txt or tsv file")
        orig_clinvar_df = pd.read_csv(clinvar_csv_path, sep="\t")
        desired_cols = ['Name', 'GRCh38Chromosome', 'GRCh38Location']
        if not all([col in orig_clinvar_df.columns for col in desired_cols]):
            raise ValueError("Clinvar file must contain columns: Name, GRCh38Chromosome, GRCh38Location")

        # filter for variants in the transcript
        clinvar_df = orig_clinvar_df[orig_clinvar_df['Name'].str.contains(transcript_name, case=False)]
        
        def filter_out_non_missense(df, protein_change_col, keep_synonymous=False):
            # pull out protein change and drop rows without protein change (unless it's synonymous)
            # TODO FIX TO KEEP SYN
            df['protein_change'] = df[protein_change_col].str.extract(r'p\.([A-Z][a-z]{2}\d+[A-Z][a-z]{2})')
            df = df.dropna(subset=['protein_change'])

            if not keep_synonymous:
                df = df[~df['protein_change'].str.contains("=")]

            # drop non-missense (nonsense, start/stop loss)
            df = df[~df['protein_change'].str.contains("Ter")]
            df = df[~df['protein_change'].str.contains(r'Met1[A-Z][a-z]{2}')]
            return df
        
        def process_reading_frame(df, coding_change_col):
            df['coding_pos'] = df[coding_change_col].str.extract(r'(\d+)[ACTG]>[ACTG](?:\s|$)')
            df['read_frame_pos'] = (((df['coding_pos'].astype(int) - 1) % 3) +1) 
            return df

        # error out if less than three variants
        clinvar_df = filter_out_non_missense(clinvar_df, 'Name', keep_synonymous=False)
        if len(clinvar_df) < 3:
            raise ValueError("Clinvar file must have at least three variants")
        
        # create standard SNV name
        clinvar_df['snv'] = clinvar_df['Name'].str.extract(r'\d+([ACTG]>[ACTG])(?:\s|$)')
        clinvar_df = clinvar_df.dropna(subset=['snv', 'GRCh38Chromosome', 'GRCh38Location']) # keep SNV only
        clinvar_df['chr'] = clinvar_df['GRCh38Chromosome'].astype(int)
        clinvar_df['pos'] = clinvar_df['GRCh38Location'].astype(int)
        clinvar_df['ref'] = clinvar_df['snv'].apply(lambda x: x[0])
        clinvar_df['alt'] = clinvar_df['snv'].apply(lambda x: x[-1])
        clinvar_df = process_reading_frame(clinvar_df, 'Name')

        # capture gene orientation
        clinvar_df.sort_values("pos")
        if clinvar_df['coding_pos'].iloc[0] < clinvar_df['coding_pos'].iloc[-1]:
            coding_strand = "+"
        elif clinvar_df['coding_pos'].iloc[-1] < clinvar_df['coding_pos'].iloc[0]:
            coding_strand = "-"
        else:
            print("Error: gene strand could not be identified")
            raise ValueError("Alert Developer - gene strand could not be identified")
        
        if coding_strand == "-":
            # necessary as clinvar reports ref and alt on the coding strand. This code reports all variation on the + strand
            clinvar_df['ref'] = clinvar_df['ref'].apply(lambda x: str(Seq(x).reverse_complement()))
            clinvar_df['alt'] = clinvar_df['alt'].apply(lambda x: str(Seq(x).reverse_complement()))

        clinvar_df['var_id'] = clinvar_df.apply(lambda x: f"{x['chr']}_{x['pos']}_{x['ref']}_{x['alt']}", axis=1)
        # clean up df
        clinvar_df = clinvar_df[['var_id', 'chr', 'pos', 'ref', 'alt', 'protein_change', 'coding_pos', 'read_frame_pos', 'Germline classification']]
        if gnomad_csv_path:
            if not gnomad_csv_path.endswith(".csv"):
                raise ValueError("Gnomad file must be a csv file")
            gnomad_df = pd.read_csv(gnomad_csv_path)
            desired_cols = ['Allele Count', 'Allele Number', 'Allele Frequency', 'gnomAD ID']
            if not all([col in gnomad_df.columns for col in desired_cols]):
                raise ValueError("Gnomad file must contain columns: Allele Count, Allele Number, Allele Frequency, gnomAD ID")
            
            gnomad_df['var_id'] = gnomad_df['gnomAD ID'].apply(lambda x: x.replace("-", "_"))
            gnomad_for_clinvar = gnomad_df[['var_id', 'Allele Count', 'Allele Number', 'Allele Frequency']]
            clinvar_w_gnomad = pd.merge(clinvar_df, gnomad_for_clinvar, on='var_id', how='left')

            # make non-clinvar gnomad df
            gnomad_df_no_clinvar = gnomad_df[~gnomad_df['var_id'].isin(clinvar_df['var_id'])]

            # pull out protein change and drop rows without protein change
            # drop non-missense (nonsense, start/stop loss)
            gnomad_df_no_clinvar = filter_out_non_missense(gnomad_df_no_clinvar, 'Protein Consequence', keep_synonymous=True)

            gnomad_df_no_clinvar['chr'] = gnomad_df_no_clinvar['var_id'].apply(lambda x: int(x.split("_")[0]))
            gnomad_df_no_clinvar['pos'] = gnomad_df_no_clinvar['var_id'].apply(lambda x: int(x.split("_")[1]))
            gnomad_df_no_clinvar['ref'] = gnomad_df_no_clinvar['var_id'].apply(lambda x: x.split("_")[2])
            gnomad_df_no_clinvar['alt'] = gnomad_df_no_clinvar['var_id'].apply(lambda x: x.split("_")[3])
            gnomad_df_no_clinvar = process_reading_frame(gnomad_df_no_clinvar, 'Transcript Consequence')
            gnomad_df_no_clinvar = gnomad_df_no_clinvar[['var_id', 'chr', 'pos', 'ref', 'alt', 'protein_change', 'coding_pos', 'read_frame_pos','Allele Count', 'Allele Number', 'Allele Frequency']]
            var_df = pd.concat([clinvar_w_gnomad, gnomad_df_no_clinvar])
        else:
            var_df = clinvar_df

        # print df of variants not in var_df but in orig_clinvar_df (excluded variants from upload)
        #excluded_variants = orig_clinvar_df[~orig_clinvar_df['Name'].isin(var_df['Name'])]

        return var_df, coding_strand
    

    def build_desired_dfs(self):
        self.pathogenic_vars = self.var_df[self.var_df["Germline classification"].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"])]
        self.benign_vars = self.var_df[self.var_df['Germline classification'].isin(["Benign", "Likely benign", "Benign/Likely benign"])]
        self.VUS_vars = self.var_df[self.var_df['Germline classification'].isin(["Uncertain significance", "Conflicting classifications of pathogenicity"])]
        if 'Allele Count' in self.var_df:
            self.gnomad_vars = self.var_df[self.var_df['Allele Count'] >= self.allele_count_min]
        else:
            self.gnomad_vars = pd.DataFrame()
        # print the number of pathogenic, benign, and VUS variants
        # print(f"Total variants: {self.var_df.shape[0]}")
        # print(f"Pathogenic variants: {self.pathogenic_vars.shape[0]}")
        # print(f"Benign variants: {self.benign_vars.shape[0]}")
        # print(f"VUS variants: {self.VUS_vars.shape[0]}")
        # print(f"Gnomad variants: {self.gnomad_vars.shape[0]}")

        #self.unclassified_df = self.var_df[~self.var_df['var_id'].isin(self.pathogenic_vars['var_id']) & ~self.var_df['var_id'].isin(self.benign_vars['var_id']) & ~self.var_df['var_id'].isin(self.VUS_vars['var_id'])]

    def find_pam_sites(self, start, stop):
        # find all GG or CC sites in the desired region
        pam_sites = []
        for i in range(start, stop-2):
            potential_pam = self.ref_fasta[i-1:i+2] #convert to 0 indexing
            if potential_pam[1:] == "GG":
                pam_sites.append({"pam_start_loc": i, "strand": "+"})
            if potential_pam[:-1] == "CC":
                pam_sites.append({"pam_start_loc": i+2, "strand": "-"})
        return pd.DataFrame(pam_sites)
    

    def find_variant_dense_windows(self, desired_var_df):
        #TODO: in future, can search only around exons of desired gene to speed up window finding
        windows = []
        pam_sites = self.find_pam_sites(self.vars_start-20-self.rtt_len, self.vars_end+20+self.rtt_len) # add some padding to allow for spacers outside of the variant region
        for _, row in pam_sites.iterrows():
            if row['strand'] == "+":
                rtt_start = row['pam_start_loc'] - 3
                rtt_end = rtt_start + self.rtt_len - 1
            elif row['strand'] == "-":
                rtt_end = row['pam_start_loc'] + 3
                rtt_start = rtt_end - self.rtt_len+1 
            else:
                print("Error: strand not recognized. Contact developer")
                continue

            # get all variants in the window
            count = desired_var_df["pos"].between(rtt_start, rtt_end).sum()
            # count the number of variants in the df
            windows.append({'rtt_start':rtt_start, 'rtt_end':rtt_end, "num_vars":count, "peg_strand":row['strand']}) # store in 1 indexing

        # sort the windows by priority
        windows = sorted(windows, key=lambda x: x['num_vars'], reverse=True)

        # get the top, non-overlapping windows
        top_windows = []
        for window in windows:
            if not any([window['rtt_start'] <=  tw['rtt_start'] <= window['rtt_end'] or window['rtt_start'] <= tw['rtt_end'] <= window['rtt_end'] for tw in top_windows]):
                top_windows.append(window)
                #print(f"Window: {window['rtt_start']}-{window['rtt_end']}, # variants: {window['num_vars']}, strand: {window['peg_strand']}")
            if len(top_windows) == self.num_windows:
                break
        
        return top_windows


    def design_pegrnas(self, window_start, window_end, peg_strand, window_df, stop_codons=False, guide_screening_mode=False):
        # convert back to 0 indexing
        if peg_strand == "+":
            spacer_start = window_start-18 
            spacer_end = window_start+1
            spacer = self.ref_fasta[spacer_start:spacer_end+1]
            pam = self.ref_fasta[spacer_end+1:spacer_end+4]

            rtt_rev_temp = self.ref_fasta[window_start-1:window_end]
            pbs_rev = self.ref_fasta[window_start-self.pbs_len-1:window_start-1]
        elif peg_strand == "-":
            spacer_start = window_end-3
            spacer_end = window_end+16
            spacer = str(Seq(self.ref_fasta[spacer_start:spacer_end+1]).reverse_complement())
            pam = str(Seq(self.ref_fasta[spacer_start-3:spacer_start]).reverse_complement())

            rtt_rev_temp = self.ref_fasta[window_start-1:window_end]
            pbs_rev = str(Seq(self.ref_fasta[window_end:window_end+self.pbs_len]).reverse_complement())
        else:
            return 
        
        if guide_screening_mode:
            # duplicate last entry in window_df 
            no_edit_entry = window_df.iloc[[-1]]
            no_edit_entry['var_id'] = "syn_edit"
            no_edit_entry['alt'] = no_edit_entry['ref']
            window_df = no_edit_entry


        window_df["ref_seq"], window_df["alt_seq"] = zip(*window_df.apply(self.create_edit_fastas, axis=1))
        pegrna_data = {}
        for _, row in window_df.iterrows():
            local_edit_pos = row['pos'] - window_start
            if row['ref'].upper() != rtt_rev_temp[local_edit_pos]:
                print(f"Warning: Reference genome does not match the reference allele: {row['var_id']}")
                raise ValueError("WARNING ALERT DEVELOPER: ClinVar reference doesn't match hg38 fasta")
            
            if guide_screening_mode:
                rtt_rev = rtt_rev_temp
                lower_pam_changes = True
            else:
                rtt_rev = rtt_rev_temp[:local_edit_pos] + row['alt'].lower() + rtt_rev_temp[local_edit_pos+1:]
                lower_pam_changes = False

                
            if peg_strand == "-":
                rtt_rev = str(Seq(rtt_rev).reverse_complement())
            
            pegrna_data[row['var_id']] = {'strand':peg_strand, 'spacer': spacer, 'pam':pam, 'rtt':rtt_rev, 'pbs':str(Seq(pbs_rev).reverse_complement())}

        peg_df = pd.DataFrame(pegrna_data).transpose().reset_index()
        merged_df = pd.merge(window_df, peg_df, left_on='var_id', right_on='index', how='left')

        if self.disrupt_pam or guide_screening_mode:
            merged_df["rtt"], merged_df["pam_status"], merged_df["seed_status"], merged_df["aa_change"], merged_df["warnings"] = zip(*merged_df.apply(lambda x: self.disrupt_pam_seed(x["ref_seq"], x["alt_seq"], x["spacer"], x["rtt"], lower_changes=lower_pam_changes), axis=1))

        merged_df = merged_df.sort_values("pos")

        if stop_codons:
            merged_df = self.introduce_ptcs(merged_df, rtt_rev_temp, window_start, window_end, peg_strand)

        # add nicking guides
        merged_df['nicking sgrna'], merged_df['distance_to_nick'] = self.pick_nicking_guide(window_start, peg_strand)

        if guide_screening_mode:
            merged_df = merged_df.drop(columns=["chr","pos","ref","alt","protein_change","coding_pos", "Germline classification"])
            if "Allele Count" in merged_df:
                merged_df = merged_df.drop(columns=["Allele Count", "Allele Number", "Allele Frequency"])

        return merged_df


    def create_edit_fastas(self, row):
        # get the location of the variant
        pos = row["pos"]
        ref = row["ref"]
        alt = row["alt"]
        reading_frame = int(row["read_frame_pos"]) # if str(row["read_frame_pos"]) != "nan" else 2
        # get the sequence of the variant + ~60 bp on either side
        padding_input = 102 
        # ensures reading frame is maintained
        if self.coding_strand == "+":
            l_padding = padding_input + reading_frame - 1
            r_padding = padding_input - reading_frame + 3
        else:
            l_padding = padding_input - reading_frame + 3
            r_padding = padding_input + reading_frame - 1

        ref_seq = self.ref_fasta[pos - l_padding-1 : pos + r_padding]
        # add the edit to the sequence
        alt_seq = ref_seq[:l_padding] + alt + ref_seq[l_padding + 1 :]
        # create warning if ref doesn't match the reference genome
        if ref_seq[l_padding] != ref:
            print(f"Warning: Reference genome does not match the reference allele: {row['var_id']}")
            raise ValueError(f"Warning: Reference genome does not match the reference allele: {row['var_id']}")
        if (len(alt_seq)) % 3 != 0:
            print(f"Warning: Reading frame is not maintained: {row['var_id'], len(alt_seq)}")
            raise ValueError(f"Warning: Reading frame is not maintained: {row['var_id'], len(alt_seq)}")
        return ref_seq, alt_seq
    
    # TODO: DO I NEED THIS
    def find_spacer(self, ref_seq, spacer):
        # input validation: find spacer in the input
        orientation = "+"
        spacer_start = ref_seq.find(spacer)
        if spacer_start == -1:
            ref_seq = str(Seq(ref_seq).reverse_complement())
            spacer_start = ref_seq.find(spacer)
            if spacer_start == -1:
                return "spacer not found", -1, -1
            else:
                orientation = "-"
        spacer_end = spacer_start + len(spacer)

        return spacer_start, spacer_end, orientation


    def find_aa_changes(self, ref_seq, alt_seq, codon_flip):
        if len(ref_seq) % 3 != 0:
            ref_seq = ref_seq[:-(len(ref_seq) % 3)]
        if len(alt_seq) % 3 != 0:
            alt_seq = alt_seq[:-(len(alt_seq) % 3)]
        # split into codons
        ref_codons = [ref_seq[i : i + 3].upper() for i in range(0, len(ref_seq), 3)]
        alt_codons = [alt_seq[i : i + 3].upper() for i in range(0, len(alt_seq), 3)]

        if codon_flip:
            ref_codons = [str(Seq(codon).reverse_complement()) for codon in ref_codons]
            alt_codons = [str(Seq(codon).reverse_complement()) for codon in alt_codons]

        # count the number of codons that are different
        aa_changes = []
        for x, ref_codon in enumerate(ref_codons):
            if CODON_DICT[ref_codon][0] != CODON_DICT[alt_codons[x]][0]:
                aa_changes.append(CODON_DICT[ref_codon][0] + "->"+ CODON_DICT[alt_codons[x]][0])

        return aa_changes    

    # take a reverse transcription template, find the PAM and disrupt it
    def disrupt_pam_seed(self, ref_seq, alt_seq, spacer, rtt, lower_changes=False):
        warnings = []
        seed_status = "-"
        pam_status = "-"

        # find spacer in the input
        spacer_start, spacer_end, peg_strand = self.find_spacer(ref_seq, spacer)
        
        # boolean if codons need to be reverse comp'd based on coding strand and spacer orientation
        codon_flip = ((peg_strand == "-" and self.coding_strand == "+") or (peg_strand == "+" and self.coding_strand == "-"))
        
        if peg_strand == "-":
            ref_seq = str(Seq(ref_seq).reverse_complement())
            alt_seq = str(Seq(alt_seq).reverse_complement())

        if spacer_start == "spacer not found":
            warnings.append("!!spacer not found in reference seq!!")
            return rtt, pam_status, seed_status, str(self.find_aa_changes(ref_seq, alt_seq, codon_flip)), warnings

        # check if pam has been disrupted by edit
        ref_pam = ref_seq[spacer_end:spacer_end + 3]
        alt_pam = alt_seq[spacer_end:spacer_end + 3]
        if ref_pam[1:] != "GG":
            warnings.append("PAM not found in reference seq: issue with spacer")
            return rtt, pam_status, seed_status, str(self.find_aa_changes(ref_seq, alt_seq, codon_flip)), warnings
        if alt_pam[1:] != "GG":
            pam_status = "edit disrupts the PAM: " + alt_pam
            return rtt, pam_status, seed_status, str(self.find_aa_changes(ref_seq, alt_seq, codon_flip)), warnings

        # check if seed has been disrupted by edit
        ref_seed = ref_seq[spacer_end - 3 : spacer_end]
        alt_seed = alt_seq[spacer_end - 3 : spacer_end]
        if ref_seed != alt_seed:
            seed_status = "edit disrupts the seed: " + ref_seed + "->" + alt_seed
            return rtt, pam_status, seed_status, str(self.find_aa_changes(ref_seq, alt_seq, codon_flip)), warnings

        spacer_pam = spacer + ref_pam
        pam_end = spacer_start + len(spacer_pam)

        # split the new gdna into codons
        codons = [alt_seq[i : i + 3] for i in range(0, len(alt_seq), 3)]
        if codon_flip:
            codons = [str(Seq(codon).reverse_complement()) for codon in codons]

        # find the codons that the PAM and seed are in
        reading_frame = pam_end % 3
        if reading_frame == 0:
            pam_codons = [codons[(pam_end // 3) - 1]]
            seed_codons = [codons[(pam_end // 3) - 2]]
        else:
            pam_codons = [codons[(pam_end // 3) - 1], codons[pam_end // 3]]
            seed_codons = [codons[(pam_end // 3) - 2], codons[(pam_end // 3) - 1]]
        
        # make all synonymous changes to the PAM and seed codons
        syn_codon_changes = make_synonymous_changes(pam_codons + seed_codons)

        # find each synonymous change that disrupts the PAM
        possible_pam_disrupts = []

        for x, orig_codon in enumerate(pam_codons):
            for change in syn_codon_changes[orig_codon]:
                # flip back to original orientation if previously flipped
                syn_codon = str(Seq(change["syn_codon"]).reverse_complement()) if codon_flip else change["syn_codon"]
                if len(pam_codons) == 1:
                    new_pam_region = syn_codon
                    old_pam_codon_region = orig_codon
                else:
                    pam_codon1 = str(Seq(pam_codons[0]).reverse_complement()) if codon_flip else pam_codons[0]
                    pam_codon2 = str(Seq(pam_codons[1]).reverse_complement()) if codon_flip else pam_codons[1]
                    new_pam_region = syn_codon + pam_codon2 if x == 0 else pam_codon1 + syn_codon
                    old_pam_codon_region = pam_codon1 + pam_codon2
                
                pam_start_index = reading_frame
                new_pam = new_pam_region[pam_start_index:pam_start_index+3]

                if new_pam[1:] != "GG":
                    possible_pam_disrupts.append(
                        {
                            "old_pam_codon_region": old_pam_codon_region,
                            "new_pam_codon_region": new_pam_region,
                            "ref_codon": orig_codon,
                            "alt_codon": change["syn_codon"],
                            "old_pam": spacer_pam[-3:],
                            "new_pam": new_pam,
                            "alt_codon_freq": change["freq"],
                            "freq_diff": abs(change["freq"] - CODON_DICT[orig_codon][2]),
                        }
                    )

        # find each synonymous change that disrupts the seed by at least one
        # limitation: can only disrupt one codon right now
        possible_seed_disrupts = []
        for x, orig_codon in enumerate(seed_codons):
            for change in syn_codon_changes[orig_codon]:
                syn_codon = str(Seq(change["syn_codon"]).reverse_complement()) if codon_flip else change["syn_codon"]
                if len(seed_codons) == 1:
                    new_seed_region = syn_codon
                    old_seed_codon_region = orig_codon
                else:
                    seed_codon1 = str(Seq(seed_codons[0]).reverse_complement()) if codon_flip else seed_codons[0]
                    seed_codon2 = str(Seq(seed_codons[1]).reverse_complement()) if codon_flip else seed_codons[1]
                    new_seed_region = syn_codon + seed_codon2 if x == 0 else seed_codon1 + syn_codon
                    old_seed_codon_region = seed_codon1 + seed_codon2

                seed_start_index = reading_frame
                new_seed = new_seed_region[seed_start_index:seed_start_index+3]

                # if any edits made before the seed, ignore this change (out of prime editing window)
                if count_mismatch(new_seed_region[:seed_start_index], old_seed_codon_region[:seed_start_index]) > 0:
                    continue

                num_seed_mismatch = count_mismatch(new_seed, ref_seed)
                if num_seed_mismatch >= 1:
                    possible_seed_disrupts.append(
                        {
                            "old_seed_codon_region": old_seed_codon_region[seed_start_index:],
                            "new_seed_codon_region": new_seed_region[seed_start_index:],
                            "ref_codon": orig_codon,
                            "alt_codon": change["syn_codon"],
                            "old_seed": ref_seed,
                            "new_seed": new_seed,
                            "num_mismatch": num_seed_mismatch,
                            "alt_codon_freq": change["freq"],
                            "freq_diff": abs(change["freq"] - CODON_DICT[orig_codon][2]),
                        }
                    )

        # get ready to modify pam/seed
        rtt_new = Seq(rtt)
        rtt_pam = rtt_new[3:6]
        rtt_seed = rtt_new[0:3]

        def swap_pam(possible_pam_disrupts, rtt_new):
            # choose edited synonymous codon with the highest frequency
            pam_edit = min(possible_pam_disrupts, key=lambda x: x["freq_diff"])
            if len(pam_edit["new_pam_codon_region"]) == 3:
                rtt_new = rtt_new[:3] + pam_edit["new_pam_codon_region"] + rtt_new[6:]
            else: # length 6
                # find the new pam region
                start = rtt_new[1:8].upper().find(pam_edit["old_pam_codon_region"]) + 1
                if start == 0:
                    print("Error: old pam not found in rtt -- Alert Developer")
                    raise ValueError("Error: old pam not found in rtt -- Alert Developer")
                else:
                    rtt_new = rtt_new[:start] + pam_edit["new_pam_codon_region"] + rtt_new[start+6:]

            return rtt_new, pam_edit["new_pam"]

        def swap_seed(possible_seed_disrupts, rtt_new):
            seed_edit = min(possible_seed_disrupts, key=lambda x: x["freq_diff"])
            if len(seed_edit["new_seed_codon_region"]) == 3:
                rtt_new = seed_edit["new_seed_codon_region"] + rtt_new[3:]
            else: # length 6
                # find the new pam region
                start = rtt_new[0:5].upper().find(seed_edit["old_seed_codon_region"])
                if start == -1:
                    print("Error: old seed not found in rtt -- Alert Developer")
                    raise ValueError("Error: old seed not found in rtt -- Alert Developer")
                else:
                    rtt_new = rtt_new[:start] + seed_edit["new_seed_codon_region"] + rtt_new[start+len(seed_edit["new_seed_codon_region"]):]

            return rtt_new, seed_edit["new_seed"]

        # try to disrupt the PAM
        if rtt_pam == spacer_pam[-3:]:
            if len(possible_pam_disrupts) > 0:
                rtt_new, new_pam = swap_pam(possible_pam_disrupts, rtt_new)
                pam_status = rtt_pam + " -> " + new_pam
            else:
                # no synonymous changes disrupt the PAM
                pam_status = "-"
        else:
            # edit made a new pam (already checked above that pam wasn't disrupted by edit)
            pam_status = "edit created new PAM"

        # try to disrupt the seed if necessary and get seed status
        if rtt_seed == spacer_pam[-6:-3]:
            if pam_status in ["-", "edit created new PAM"]:
                if len(possible_seed_disrupts) > 0:
                    rtt_new, new_seed = swap_seed(possible_seed_disrupts, rtt_new)
                    seed_status = rtt_seed + " -> " + new_seed
                else:
                    seed_status = "-"
                    warnings.append("No synonymous changes disrupt the seed or PAM")

        # before returning the new rtt, confirm that only one AA was changed
        new_alt_seq = ref_seq[:spacer_end - 3] + rtt_new + ref_seq[spacer_end-3 + len(rtt_new):]
        aa_changes = self.find_aa_changes(ref_seq, new_alt_seq, codon_flip)
        if len(aa_changes) >1:
            warnings.append("More than one AA change in the edited sequence")
        if len(aa_changes) == 0:
            warnings.append("No AA change in the edited sequence")
        
        if lower_changes:
            rtt_lower = ""
            for x, rtt_new_base in enumerate(rtt_new):
                rtt_lower += rtt_new_base.lower() if rtt_new_base != rtt[x] else rtt_new_base
            rtt_new = rtt_lower
        
        return str(rtt_new), pam_status, seed_status, str(aa_changes), warnings

    def find_first_codon(self, peg_df, start, end, codon_flip, strand):
        if strand == "+":
            first_rtt_variant = peg_df.iloc[0]
            first_var_location = first_rtt_variant['pos'] - start
            last_rtt_variant = peg_df.iloc[-1]
        else:
            first_rtt_variant = peg_df.iloc[-1]
            first_var_location = end - first_rtt_variant['pos']
            last_rtt_variant = peg_df.iloc[0]

        read_frame_pos = first_rtt_variant["read_frame_pos"]
        
        if codon_flip:
            read_frame_pos = {1:3, 2:2, 3:1}[read_frame_pos]
        
        first_codon_idx = first_var_location - read_frame_pos + 1
        
        while first_codon_idx < 0 or first_codon_idx > 2:
            if first_codon_idx < 0:
                first_codon_idx +=3
            if first_codon_idx > 2:
                first_codon_idx -= 3
        
        return first_codon_idx, last_rtt_variant

    def introduce_ptcs(self, merged_df,rtt_rev_temp, window_start, window_end, strand):
        # find first codon in rtt
        codon_flip = ((strand == "-" and self.coding_strand == "+") or (strand == "+" and self.coding_strand == "-"))
        first_codon_idx, last_rtt_variant = self.find_first_codon(merged_df, window_start, window_end, codon_flip, strand)

        if strand == "-":
            rtt_rev_temp = str(Seq(rtt_rev_temp).reverse_complement())

        #pam disruption rtt
        pam_disrupt_rtt = last_rtt_variant['rtt'][:8]
            
        index = first_codon_idx
        ptc_num = 1
        stop_rtts = []
        while index+2 < self.rtt_len:
            if codon_flip:
                rtt_rev = rtt_rev_temp[:index] + "tca" + rtt_rev_temp[index+3:]
            else:
                rtt_rev = rtt_rev_temp[:index] + "tga" + rtt_rev_temp[index+3:]

            if index > 7:
                rtt_rev = pam_disrupt_rtt + rtt_rev[8:]
                if 'pam_status' in last_rtt_variant:
                    pam_status = last_rtt_variant['pam_status']
                    seed_status = last_rtt_variant['seed_status']
            else:
                pam_status = "stop codon disrupts PAM/seed region"
                seed_status = "-"
                 
            pos = window_start + index if strand == "+" else window_end - index

            # check for aa changes:
            aa_changes = self.find_aa_changes(rtt_rev_temp[first_codon_idx:], rtt_rev[first_codon_idx:], codon_flip)
            warnings = ["error with ptc introduction"] if (len(aa_changes) != 1 or '->End' not in aa_changes[0]) else []
            if len(aa_changes) == 0:
                warnings.append("possible endogenous stop codon site")

            ptc_data = {'var_id': "PTC_" + str(ptc_num)+"_"+str(pos), 'chr': last_rtt_variant['chr'], 'pos': pos, 'strand':strand, 'spacer': last_rtt_variant['spacer'], 'pam':last_rtt_variant['pam'], 'rtt':rtt_rev, 'pbs':last_rtt_variant['pbs']}
            if 'pam_status' in last_rtt_variant:
                ptc_data.update({'pam_status': pam_status, 'seed_status': seed_status, 'aa_change': str(aa_changes), 'warnings': warnings})
                
            stop_rtts.append(ptc_data)
            index += 3
            ptc_num +=1
        
        stop_df = pd.DataFrame(stop_rtts)
        merged_df = pd.concat([merged_df, stop_df])
        merged_df = merged_df.sort_values("pos")

        return merged_df

    def pick_nicking_guide(self, pe_nick_site, pe_nick_strand, min_dist = 40, max_dist = 100):
        # find all GG or CC sites in the desired region
        pam_sites = self.find_pam_sites(pe_nick_site-max_dist, pe_nick_site-min_dist)
        pam_sites = pd.concat([pam_sites, self.find_pam_sites(pe_nick_site+min_dist, pe_nick_site+max_dist)])
        if pe_nick_strand == "+":
            pam_sites = pam_sites[pam_sites["strand"] == "-"]
        elif pe_nick_strand == "-":
            pam_sites = pam_sites[pam_sites["strand"] == "+"]
        else:
            return "strand not recognized"

        # find the closest site
        if pam_sites.shape[0] == 0:
            return f"no ngrna site found {min_dist}-{max_dist}bp away from pe nick site", 0
        pam_sites.reset_index(drop=True, inplace=True)
        closest_site = pam_sites.iloc[pam_sites["pam_start_loc"].sub(pe_nick_site).abs().idxmin()]
        
        if closest_site['strand'] == "+":
            sgrna = self.ref_fasta[closest_site["pam_start_loc"]-21 : closest_site["pam_start_loc"]-1]
        elif closest_site['strand'] == "-":
            sgrna = str(Seq(self.ref_fasta[closest_site["pam_start_loc"] : closest_site["pam_start_loc"]+20]).reverse_complement())
        else:
            return

        return sgrna, closest_site["pam_start_loc"] - pe_nick_site

    def clipe_prepare_oligo_df(self, final_peg_df):
        # export for IDT -- PCR strategy
        idt_df = pd.DataFrame()
        idt_df["var_id"] = final_peg_df["var_id"]
        idt_df["5' flanking left"] = "ggctac"
        idt_df["5' Bsa1 recognition site"] = "ggtctcc"
        pre_spacer_seq = final_peg_df["spacer"].apply(lambda x: "cacc" if x[0].lower()=="g" else "caccg")
        idt_df["idt_spacer"] = pre_spacer_seq + final_peg_df["spacer"] + "gtttt"
        idt_df["scaffold"] = "AGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG".lower()
        idt_df["idt_rtt"] = final_peg_df["rtt"].apply(lambda x: "gtgc" + str(Seq(x).reverse_complement()))
        idt_df['idt_pbs'] = final_peg_df['pbs']
        idt_df["3' Bsa1 recognition site"] = "cgcgtgagacc"
        idt_df["3' flanking right"] = "gtagcc"

        idt_df["full_peg"] = (
            idt_df["5' flanking left"]
            + idt_df["5' Bsa1 recognition site"]
            + idt_df["idt_spacer"]
            + idt_df["scaffold"]
            + idt_df["idt_rtt"]
            + idt_df["idt_pbs"]
            + idt_df["3' Bsa1 recognition site"]
            + idt_df["3' flanking right"]
        )
        idt_df['PCR_F'] = idt_df["full_peg"].apply(lambda x: x[:20])
        idt_df['PCR_R'] = idt_df["full_peg"].apply(lambda x: str(Seq(x[-20:]).reverse_complement()))

        warnings = []
        # Confirm no internal BsaI sites
        for seq in idt_df["full_peg"]:
            seq = seq.upper()
            num_sites = seq.count("GGTCTC") + seq.count("GAGACC")
            if num_sites != 2:
                warnings.append(f"{num_sites} BsaI sites: likely internal BsaI site")
            else:
                warnings.append("")
        idt_df["oligo_warnings"] = warnings

        return idt_df

def build_files_for_jellyfish(final_df, screening_df=None):
    # take the rtt and var_ids from the final_df and prep for fasta output
    rtt_df = final_df[["var_id", "rtt"]]
    if screening_df is not None:
        rtt_df = pd.concat([rtt_df, screening_df[["var_id", "rtt"]]])
    rtt_df = rtt_df.sort_values("var_id")
    # check for duplicate var or rtts
    if rtt_df["var_id"].duplicated().any():
        print("Duplicate var_ids in the final_df")
        raise ValueError("Duplicate var_ids in the final_df")
    if rtt_df["rtt"].duplicated().any():
        print("Duplicate rtts in the final_df")
        raise ValueError("Duplicate rtts in the final_df")
    
    # write to txt variable that can later be written to a file
    fasta_str = ""
    for _, row in rtt_df.iterrows():
        fasta_str += f">{row['var_id']}\n{row['rtt']}\n"

    return rtt_df, fasta_str
    
def prep_anneal_df(file_prefix, oligo_df, mode="nick", overhangs=['cacc', 'aaac']):
    oligo_df['name'] = oligo_df['editing_window'].apply(lambda x: f"{file_prefix}window_{x}_{mode}")
    if mode in ["nick", "screen_spacer"]:
        oligo_df['oligo'] = oligo_df['oligo'].apply(lambda x: "g" + x if x[0].upper() != "G" else  x)
    
    rev_comp_df = oligo_df.copy()
    rev_comp_df['oligo'] = rev_comp_df['oligo'].apply(lambda x: str(Seq(x).reverse_complement()))
    oligo_df['name'] = oligo_df['name'].apply(lambda x: f"{x}_top")
    oligo_df['oligo'] = oligo_df['oligo'].apply(lambda x: overhangs[0] + x)
    rev_comp_df['name'] = rev_comp_df['name'].apply(lambda x: f"{x}_bottom")
    rev_comp_df['oligo'] = rev_comp_df['oligo'].apply(lambda x: overhangs[1] + x)

    full_oligo_df = pd.concat([oligo_df, rev_comp_df])
    full_oligo_df['scale'] = "25nm"
    full_oligo_df['purification'] = "STD"
    full_oligo_df.sort_values("editing_window", inplace=True)

    return full_oligo_df

def prep_nicking_order_df(file_prefix, peg_df):
    nick_df = peg_df[['nicking sgrna', 'editing_window']].drop_duplicates()
    nick_df.columns = ['oligo', 'editing_window']
    nick_df = prep_anneal_df(file_prefix, nick_df, mode="nick", overhangs=['cacc', 'aaac'])
    return nick_df

def prep_screening_order_df(file_prefix, screening_df):
    overhangs = {'spacer': ['cacc', 'ctct'], 'extension': ['gtgc', 'cgcg']}
    screening_df = screening_df[['spacer', 'rtt', 'pbs', 'editing_window']]
    screening_df['extension'] = screening_df.apply(lambda x: str(Seq(x.rtt).reverse_complement()).upper() + x.pbs, axis=1)
    screen_oligo_df = pd.DataFrame()
    for col in ['spacer', 'extension']:
        sub_df = screening_df[[col, 'editing_window']]
        sub_df.columns = ['oligo', 'editing_window']
        if col == 'spacer':
            sub_df['oligo'] = sub_df['oligo'].apply(lambda x: x + "gtttt")
        temp_screening_df = prep_anneal_df(file_prefix, sub_df, mode=f"screen_{col}", overhangs=overhangs[col])
        screen_oligo_df = pd.concat([screen_oligo_df, temp_screening_df])

    return screen_oligo_df


CODON_DICT = {
        'GGG':['Gly','G', 0.25],'GGA':['Gly','G', 0.25],'GGT':['Gly','G', 0.16],'GGC':['Gly','G', 0.34],
        'GAG':['Glu','E', 0.58],'GAA':['Glu','E', 0.42],'GAT':['Asp','D', 0.46],'GAC':['Asp','D', 0.54],
        'GTG':['Val','V', 0.47],'GTA':['Val','V', 0.11],'GTT':['Val','V', 0.18],'GTC':['Val','V', 0.24],
        'GCG':['Ala','A', 0.11],'GCA':['Ala','A', 0.23],'GCT':['Ala','A', 0.26],'GCC':['Ala','A', 0.4],
        'AGG':['Arg','R', 0.2],'AGA':['Arg','R', 0.2],'AGT':['Ser','S', 0.15],'AGC':['Ser','S', 0.24],
        'AAG':['Lys','K', 0.58],'AAA':['Lys','K', 0.42],'AAT':['Asn','N', 0.46],'AAC':['Asn','N', 0.54],
        'ATG':['Met','M', 1],'ATA':['Ile','I', 0.16],'ATT':['Ile','I', 0.36],'ATC':['Ile','I', 0.48],
        'ACG':['Thr','T', 0.12],'ACA':['Thr','T', 0.28],'ACT':['Thr','T', 0.24],'ACC':['Thr','T', 0.36],
        'TGG':['Trp','W', 1],'TGA':['End','X', 0.52],'TGT':['Cys','C', 0.45],'TGC':['Cys','C', 0.55],
        'TAG':['End','X', 0.2],'TAA':['End','X', 0.28],'TAT':['Tyr','Y', 0.43],'TAC':['Tyr','Y', 0.57],
        'TTG':['Leu','L', 0.13],'TTA':['Leu','L', 0.07],'TTT':['Phe','F', 0.45],'TTC':['Phe','F', 0.55],
        'TCG':['Ser','S', 0.06],'TCA':['Ser','S', 0.15],'TCT':['Ser','S', 0.18],'TCC':['Ser','S', 0.22],
        'CGG':['Arg','R', 0.21],'CGA':['Arg','R', 0.11],'CGT':['Arg','R', 0.08],'CGC':['Arg','R', 0.19],
        'CAG':['Gln','Q', 0.75],'CAA':['Gln','Q', 0.25],'CAT':['His','H', 0.41],'CAC':['His','H', 0.59],
        'CTG':['Leu','L', 0.41],'CTA':['Leu','L', 0.07],'CTT':['Leu','L', 0.13],'CTC':['Leu','L', 0.2],
        'CCG':['Pro','P', 0.11],'CCA':['Pro','P', 0.27],'CCT':['Pro','P', 0.28],'CCC':['Pro','P', 0.33],
    }

def count_mismatch(seq1, seq2):
        return sum([1 for x in range(len(seq1)) if seq1[x] != seq2[x]])

def make_synonymous_changes(codons):
    codon_changes = {}
    for codon in codons:
        codon_changes[codon] = []
        for key in CODON_DICT:
            if CODON_DICT[key][0] == CODON_DICT[codon][0] and key != codon:
                codon_changes[codon].append({'syn_codon': key, 'aa': CODON_DICT[key][0], 'freq': CODON_DICT[key][2]})
    return codon_changes
