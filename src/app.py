import faicons as fa
from pathlib import Path
import pandas as pd
import plotly.graph_objs as go
from datetime import date
from clipe_guide_design_methods import *

# Load data and compute static values
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_plotly

app_dir = Path(__file__).parent
gene_data = pd.read_csv(app_dir / "genome_files/hg38_transcripts.tsv", sep="\t")
gene_data['transcript_id'] = gene_data['transcript_id'].apply(eval)
gene_names = list(gene_data["gene_id"])
gene_data_dict = gene_data.set_index("gene_id")['transcript_id'].to_dict()


ICONS = {
    "windows": fa.icon_svg("window-restore", "regular"),
    "variants": fa.icon_svg("dna"),
    "vus": fa.icon_svg("circle-question"),
    "plp": fa.icon_svg("circle-check"),
    "blb": fa.icon_svg("circle-check"),
    "gnomad": fa.icon_svg("people-group"),
    "ellipsis": fa.icon_svg("ellipsis"),
}

# Add page title and sidebar
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_radio_buttons(
         "design_strategy",
         "Design strategy:",
         {
             "vus": "Introduce missense VUS variants",
             "pos_control": "Introduce missense PLP/BLB variants",
         },
         width="100%",
        ),
        #ui.input_text("GeneName", label="Gene Name", value="TSC2"),
        ui.input_selectize("gene", label="Gene Name", choices=[], width="100%"),
        #ui.input_text("transcript", label="Transcript (Refseq)", value="NM_000548.5"),
        ui.input_select("transcript", label="Transcript (Refseq)", choices=[], width="100%"),
        ui.input_numeric("num_designs", label="Number of epegRNA libraries to design", value=12, step=1),
        ui.input_numeric("length_rtt", label="RTT Length", value=40, step=1),
        ui.input_numeric("length_pbs", label="PBS Length", value=10, step=1),
        ui.input_radio_buttons("disrupt_pam", "Disrupt PAM/Seed region with synonymous variants?", {"yes": "Yes", "no": "No"}, width="100%"),
        ui.input_file("clinvar_csv", "Choose missense ClinVar csv to upload:", multiple=False),
        ui.input_file("gnomad_csv", "Choose gnomAD csv to upload:", multiple=False),
        ui.input_numeric("allele_min", label="Minimum allele count for gnomad (if applicable)", value=5, step=1),
        ui.input_checkbox_group(  
            "checkbox_group",  
            "Include/exclude certain variant types",  
            {  
                "BLB": "Truth_Set_BLB (Clinvar)",  
                "PLP": "Truth_Set_PLP (Clinvar)",  
                "GNOMAD": "gnomAD variants (synonymous and missense; gnomAD)",
                "PTC": "PTC variants to induce LoF for assay validation/calibration",  
            },
            width="100%"),
        ui.input_action_button("action_button", "Design epegRNA Libraries"),
        open="desktop",
        width=400
    ),
    ui.layout_columns(
        ui.value_box(
            "Editing Windows", ui.output_ui("total_windows"), showcase=ICONS["windows"], theme = "primary"
        ),
        ui.value_box(
            "pegRNA designs", ui.output_ui("total_pegs"), showcase=ICONS["variants"], theme = "primary", 
        ),
        ui.value_box(
            "VUS variants", ui.output_ui("total_vus"), showcase=ICONS["vus"], theme = "primary"
        ),
        ui.value_box(
            "PLP variants", ui.output_ui("total_plp"), showcase=ICONS["plp"], theme = "primary"
        ),
        ui.value_box(
            "BLB variants", ui.output_ui("total_blb"), showcase=ICONS["blb"], theme = "primary"
        ),
        ui.value_box(
            "Gnomad variants", ui.output_ui("total_gnomad"), showcase=ICONS["gnomad"], theme = "primary", 
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
        #TODO INSERT LINKS
        id="download_area",
        ),
        ui.card(
            ui.card_header(
                "pegRNA distribution",
                class_="d-flex justify-content-between align-items-center",
            ),
            output_widget("peg_dist_chart"),
            #full_screen=True,
        ),
        col_widths=[6, 6, 12],
    ),
    ui.include_css(app_dir / "styles.css"),
    title="CliPE pegRNA Designer",
    fillable=True,
)


def server(input, output, session):
    ui.update_selectize("gene", choices=gene_names, selected="TSC2", server=True)

    @reactive.effect
    @reactive.event(input.gene)
    def _():
        if input.gene() == "":
            choices = []
        else:
            choices = gene_data_dict[input.gene()]
            ui.update_select("transcript", choices=choices)
    
    @reactive.effect
    @reactive.event(input.example)
    def _():
        #TODO: example input with default csvs
        ui.update_checkbox_group("checkbox_group", selected=["BLB", "PLP", "PTC"])
    
    

    @render.data_frame
    @reactive.event(input.action_button)
    def build_peg_df():
        # remove download button if it exists
        ui.remove_ui("#download")
        global peg_df
        peg_df = pd.DataFrame()
        # validate input files exist
        if not input.clinvar_csv():
            raise ValueError("Clinvar file must be provided")
        if not input.gnomad_csv():
            gnomad_file = None
        else:
            gnomad_file = input.gnomad_csv()[0]["datapath"]
        
        if input.disrupt_pam() == "yes":
            disrupt_pam = True
        else:
            disrupt_pam = False
    
        expt = clipe_expt(input.transcript().split(" ")[0], input.clinvar_csv()[0]["datapath"], gnomad_file, input.length_pbs(), input.length_rtt(), input.num_designs(), disrupt_pam, input.design_strategy(), input.checkbox_group(), input.allele_min())
        peg_df = expt.run_guide_design()

        ui.insert_ui(
            ui.download_button("download", "Download CSV"),
            selector="#download_area",
            where="beforeEnd",
        )

        return peg_df

    @render.download()
    def download():
        path = Path(__file__).parent / f"{input.gene()}_clipe_peg_df_{date.today()}.csv"
        print(path)
        peg_df.to_csv(path, index=False)
        return str(path)


    @render.ui
    @reactive.event(input.action_button)
    def total_windows():
        if peg_df.shape[0] > 0:
            windows = peg_df["editing_window"].nunique()
            return windows

    @render.ui
    @reactive.event(input.action_button)
    def total_pegs():
        if peg_df.shape[0] > 0:
            return peg_df.shape[0]

    @render.ui
    @reactive.event(input.action_button)
    def total_vus():
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Uncertain significance", "Conflicting classifications of pathogenicity"])].shape[0]
        
    @render.ui
    @reactive.event(input.action_button)
    def total_plp():
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"])].shape[0]
        
    
    @render.ui
    @reactive.event(input.action_button)
    def total_blb():
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Benign", "Likely benign", "Benign/Likely benign"])].shape[0]

    @render.ui
    @reactive.event(input.action_button)
    def total_gnomad():
        if peg_df.shape[0] > 0:
            if 'Allele Count' in peg_df:
                return peg_df[peg_df['Allele Count'] >= input.allele_min()].shape[0]
            else:
                return 0

    @render_plotly
    @reactive.event(input.action_button)
    def peg_dist_chart():
    # Generate a random signal
        x = list(peg_df["coding_pos"])
        x = [int(i) for i in x]
        layout = go.Layout(
            yaxis = dict(range=[-.5, .5], showticklabels=False),  # Set y-axis scale from 0 to 1
            xaxis = dict(range=[0, max(x)+10]),
            plot_bgcolor='white',  # Change background to white
            title=input.gene() + " PE windows",
        )
        # Plot the chart
        fig = go.Figure(layout=layout)

        fig.add_shape(
            type="rect",
            x0=0,
            y0=-.3,
            x1=max(x)+10,
            y1=.3,
            line_width=3,
            opacity=0.8,
        )
        windows = peg_df.groupby("editing_window").agg(Start=("coding_pos", "min"), Finish=("coding_pos", "max")).reset_index()
        for _, window in windows.iterrows():
            fig.add_shape(
                type="rect",
                x0=window['Start'],
                y0=-0.3,
                x1=window['Finish'],
                y1=0.3,
                fillcolor="Red",
                opacity=.6,
                line_width=0,
                label=dict(text=str(window['editing_window']), font=dict(color="White"))
            )
        
        return fig

app = App(app_ui, server)