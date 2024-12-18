# import faicons as fa
from pathlib import Path
import pandas as pd
import plotly.graph_objs as go
from datetime import date
from clipe_guide_design_methods import *
import tempfile
import shutil

# Load data and compute static values
from shiny import App, reactive, render, ui, req
from shinywidgets import output_widget, render_plotly

app_dir = Path(__file__).parent
gene_data = pd.read_csv(app_dir / "genome_files/hg38_transcripts.tsv", sep="\t")
gene_data['transcript_id'] = gene_data['transcript_id'].apply(eval)
gene_names = list(gene_data["gene_id"])
gene_data_dict = gene_data.set_index("gene_id")['transcript_id'].to_dict()

cds_data = pd.read_csv(app_dir / "genome_files/hg38_transcript_coords.tsv", sep="\t")
cds_data['cds_lengths'] = cds_data['cds_lengths'].apply(eval)
cds_data['cds_exons'] = cds_data['cds_exons'].apply(eval)

example_clinvar_path = str(app_dir / "example_input/clinvar_result.txt")
example_gnomad_path = str(app_dir / "example_input/gnomAD_v4.1.0_ENSG00000103197_2024_11_03_20_28_38.csv")


# ICONS = {
#     "windows": fa.icon_svg("window-restore", "regular"),
#     "variants": fa.icon_svg("dna"),
#     "vus": fa.icon_svg("circle-question"),
#     "plp": fa.icon_svg("circle-check"),
#     "blb": fa.icon_svg("circle-check"),
#     "gnomad": fa.icon_svg("people-group"),
#     "ellipsis": fa.icon_svg("ellipsis"),
# }

# Add page title and sidebar
app_ui = ui.page_navbar(
    ui.nav_panel("pegRNA Designer", 
        ui.page_sidebar(
        ui.sidebar(
            ui.input_radio_buttons(
            "design_strategy",
            ui.span("Design strategy:", style="text-decoration: underline;"),
            {
                "vus": "Introduce missense VUS variants",
                "pos_control": "Introduce missense PLP/BLB variants",
            },
            width="100%",
            ),
            ui.input_selectize("gene", label="Gene Name", choices=[], width="100%"),
            ui.input_select("transcript", label="Transcript (Refseq)", choices=[], width="100%"),
            ui.input_numeric("num_designs", label="Number of epegRNA libraries to design", value=12, step=1, min=1),
            ui.input_file("clinvar_csv", "Choose missense ClinVar csv to upload:", multiple=False),
            ui.input_file("gnomad_csv", "Choose gnomAD csv to upload:", multiple=False),
            ui.input_checkbox_group(  
                "checkbox_group",  
                ui.span("Include additional variants in editing windows:", style="text-decoration: underline;"),  
                {  
                    "BLB": "Benign/Likely Benign variants",  
                    "PLP": "Pathogenic/Likely Pathogenic variants",  
                    "GNOMAD": "gnomAD variants with allele count >= 5",
                    "PTC": "PTC variants to induce LoF for assay validation/calibration",  
                },
                width="100%"),
            ui.span(ui.input_action_button("action_button", "Design epegRNA Libraries"), style="align-self: center;"),
            ui.accordion(ui.accordion_panel('Additional Options',
                ui.input_numeric("length_rtt", label="RTT Length", value=40, step=1, min=0, max=60),
                ui.input_numeric("length_pbs", label="PBS Length", value=10, step=1, min=4, max=20),
                ui.input_radio_buttons("disrupt_pam", "Disrupt PAM/Seed region? (synonymous variants)", {"yes": "Yes", "no": "No"}, width="100%"),
                ui.input_numeric("allele_min", label="gnomAD minimum allele count", value=5, step=1, min=0),
                ),id="additional_options", multiple=False, open=False),
            open="desktop",
            width=375,
            padding=20,
            fillable=True,
        ),
        ui.layout_columns(
            ui.value_box(
                "# pegRNA designs", ui.output_ui("total_pegs"),  theme = ui.value_box_theme(bg = "#e6f2fd", fg = "#0B538E")
            ),
            ui.value_box(
                "# VUS variants", ui.output_ui("total_vus"),  theme = ui.value_box_theme(bg = "#e6f2fd", fg = "#0B538E")
            ),
            ui.value_box(
                "# PLP variants", ui.output_ui("total_plp"),  theme= ui.value_box_theme(bg = "#e6f2fd", fg = "#0B538E")
            ),
            ui.value_box(
                "# BLB variants", ui.output_ui("total_blb"),  theme = ui.value_box_theme(bg = "#e6f2fd", fg = "#0B538E")
            ),
            fill=False,
        ),
        ui.layout_columns(
            ui.card(
                ui.card_header("pegRNA designs"), ui.output_data_frame("build_peg_df"), full_screen=True
            ),
            ui.card(
                ui.card_header(
                    "Downloads",
                    class_="d-flex justify-content-between align-items-center",
                ),
                ui.output_ui("download_area"),
            ),
            ui.card(
                ui.card_header(
                    "pegRNA distribution",
                    class_="d-flex justify-content-between align-items-center",
                ),
                output_widget("peg_dist_chart"),
            ),
            col_widths=[8, 4, 12],
            row_heights=(2,1),
            fillable=True,
        ),
        ui.include_css(app_dir / "styles.css"),
        fillable=True,
        ),
    ),
    ui.nav_spacer(),
    ui.nav_control(ui.a('CliPE Home', href='http://calhoujd.github.io', target="_blank",)),
    ui.nav_control(ui.a('Carvill Lab', href='https://sites.northwestern.edu/carvilllab/', target="_blank",)),
    ui.nav_control(ui.a(ui.img(src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png", alt="Github", height="25px"), href='https://github.com/nbodkin/clipe_pegrna_builder', target="_blank")),

    title=ui.img(src="clipe_logo.png", alt="clipe logo", height="50px"),
    window_title="CliPE pegRNA Designer",
    fillable=True,
)

def server(input, output, session):
    example_bool = reactive.value(False)
    example_count = reactive.value(0)
    peg_df_glob = reactive.value(pd.DataFrame())
    fish_df_glob = reactive.value(pd.DataFrame())
    fish_fa_txt_glob = reactive.value("")

    ui.update_selectize("gene", choices=gene_names, selected="TSC2", server=True)

    @reactive.effect
    @reactive.event(input.gene)
    def _():
        ui.remove_ui("#example_span")
        if input.gene() == "":
            choices = []
        else:
            choices = gene_data_dict[input.gene()]
            ui.update_select("transcript", choices=choices)
        
        ui.insert_ui(ui.span(ui.input_action_link("example", "Design with TSC2 example"), style="align-self: center;", id="example_span"),
                     selector="#additional_options", where="beforeBegin")
    
    @reactive.effect
    @reactive.event(input.example)
    def _():
        ui.update_selectize("gene", choices=gene_names, selected="TSC2", server=True)
        example_bool.set(True)
        example_count.set(example_count.get()+1)

    @render.data_frame
    @reactive.event(input.action_button, example_count, ignore_init=True)
    def build_peg_df():
        # remove download button if it exists
        ui.remove_ui("#download_button")
        ui.remove_ui("#download_checkbox")
        peg_df_glob.set(pd.DataFrame())
        if input.disrupt_pam() == "yes":
            disrupt_pam = True
        else:
            disrupt_pam = False
        
        # validate input files exist
        if example_bool.get():
            clinvar_file = example_clinvar_path
            gnomad_file = example_gnomad_path
        else:
            if not input.clinvar_csv():
                raise ValueError("Clinvar file must be provided")
            else:
                clinvar_file = input.clinvar_csv()[0]["datapath"]
            if not input.gnomad_csv():
                gnomad_file = None
            else:
                gnomad_file = input.gnomad_csv()[0]["datapath"]
        
        expt = clipe_expt(input.transcript().split(" ")[0], clinvar_file, gnomad_file, input.length_pbs(), input.length_rtt(), input.num_designs(), disrupt_pam, input.design_strategy(), input.checkbox_group(), input.allele_min())
        peg_df, windows = expt.run_guide_design()
        peg_df_glob.set(peg_df)

        output = expt.build_files_for_jellyfish(peg_df)
        fish_df_glob.set(output[0])
        fish_fa_txt_glob.set(output[1])
        
        # TODO: reduce sig figs in allele percentage
        ui.insert_ui(
            ui.download_button("download_button", "Download selected files", width="60%"),
            selector="#download_area",
            where="afterEnd",
        )
        ui.insert_ui(
            ui.input_checkbox_group(  
            "download_checkbox",  
            ui.span("Download the following files:", style="font-weight: bold;"),  
            {  
            "peg_tables": "pegRNA design tables (.csv)",  
            "RTTs": "RTT data for downstream analyses (.csv, .fa)",
            "idt": "IDT oPool ordering files (.xlsx)",
            },
            width="100%", selected=["peg_tables", "RTTs", "idt"]),
            selector="#download_area",
            where="afterEnd",
        )
        
        example_bool.set(False)

        return peg_df

    @render.download(filename=lambda: f"{date.today()}_{input.gene()}_clipe_designs.zip")
    def download_button():
        peg_df = peg_df_glob.get()
        fish_df = fish_df_glob.get()
        fish_fa_txt = fish_fa_txt_glob.get()
        files_to_download = input.download_checkbox()
        if len(files_to_download) > 0:
            file_prefix = f"{input.gene()}_"
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                if "peg_tables" in files_to_download:
                    peg_df.to_csv(temp_path / f"{file_prefix}full_pegRNA_designs.csv", index=False)
                if "RTTs" in files_to_download:
                    with open(temp_path / f"{file_prefix}RTT_fasta.fa", "w") as f:
                        f.write(fish_fa_txt)
                    fish_df.to_csv(temp_path / f"{file_prefix}RTT_table.csv", index=False)
                if "idt" in files_to_download:
                    idt_df = peg_df[['editing_window', "full_peg"]]
                    idt_df['editing_window'] = idt_df['editing_window'].apply(lambda x: f"{file_prefix}window_{x}")
                    idt_df.columns = ["Pool name", "Sequence"]
                    idt_df.to_excel(temp_path / f"{file_prefix}IDT_order_data.xlsx", index=False)
                zip_path = shutil.make_archive(temp_path, 'zip', temp_path)
                # yield zip file data
                with open(zip_path, "rb") as f:
                    x = f.read()
                yield(x)


    @reactive.effect
    @reactive.event(input.allele_min)
    def allele_min():
        ui.update_checkbox_group(  
            "checkbox_group",  
            choices= {  
                "BLB": "Benign/Likely Benign variants",  
                "PLP": "Pathogenic/Likely Pathogenic variants",  
                "GNOMAD": f"gnomAD variants with allele count >= {input.allele_min()}",
                "PTC": "PTC variants to induce LoF for assay validation/calibration",  
            })

    @render.ui
    def total_pegs():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] > 0:
            return peg_df.shape[0]
    
    @render.ui
    def total_vus():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Uncertain significance", "Conflicting classifications of pathogenicity"])].shape[0]
    
    @render.ui
    def total_plp():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"])].shape[0]
    
    @render.ui
    def total_blb():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Benign", "Likely benign", "Benign/Likely benign"])].shape[0]

    @render_plotly
    @reactive.event(peg_df_glob, ignore_init=True)
    def peg_dist_chart():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] != 0:
            df_to_plot = peg_df.copy().dropna(subset=["coding_pos"])
            if df_to_plot.shape[0] == 0:
                return go.Figure()
            
            x = [int(i) for i in list(df_to_plot["coding_pos"])]
            transcript_data = cds_data[cds_data['transcript_id'] == input.transcript().split(" ")[0]]
            cds_lengths = transcript_data['cds_lengths'].iloc[0]
            total_cds_length = sum(cds_lengths)
            cds_exons = transcript_data['cds_exons'].iloc[0]

            layout = go.Layout(
                yaxis = dict(range=[-.5, .5], showticklabels=False),  # Set y-axis scale from 0 to 1
                xaxis = dict(range=[1, total_cds_length]),
                plot_bgcolor='white',  # Change background to white
                title=input.gene() + " PE windows",
            )
            # Plot the chart
            fig = go.Figure(layout=layout)

            fig.add_shape(
                type="rect",
                x0=1,
                y0=-.3,
                x1=total_cds_length,
                y1=.3,
                line_width=3,
                opacity=0.8,
            )
            cds_pos = 1
            for i, cds_endpoint in enumerate(cds_lengths):
                fig.add_trace(go.Scatter(
                    x=[cds_pos, cds_pos],
                    y=[-.3, .3],
                    mode="lines",
                    line=dict(color="Gray", dash="dot"),
                    line_width=.5,
                    hoverinfo="text",
                    text=f"Exon {cds_exons[i]}",
                    showlegend=False
                ))
                cds_pos += cds_endpoint

            windows = df_to_plot.groupby("editing_window").agg(Start=("coding_pos", "min"), Finish=("coding_pos", "max")).reset_index()
            for _, window in windows.iterrows():
                fig.add_shape(
                    type="rect",
                    x0=window['Start'],
                    y0=-0.3,
                    x1=window['Finish'],
                    y1=0.3,
                    fillcolor="darkviolet",
                    opacity=.8,
                    line_width=0,
                    label=dict(text=str(window['editing_window']), font=dict(color="White"))
                )
            
            return fig

app = App(app_ui, server, static_assets=app_dir / "www")