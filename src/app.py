# import faicons as fa
from pathlib import Path
import pandas as pd
import plotly.graph_objs as go
from datetime import date
from clipe_guide_design_methods import *
import tempfile
import shutil
import gc

# Load data and compute static values
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_plotly

app_dir = Path(__file__).parent
gene_data = pd.read_csv(app_dir / "genome_files/hg38_transcripts.tsv", sep="\t")
gene_data['transcript_id'] = gene_data['transcript_id'].apply(eval)
gene_names = list(gene_data["gene_id"])
gene_data_dict = gene_data.set_index("gene_id")['transcript_id'].to_dict()

cds_data = pd.read_csv(app_dir / "genome_files/hg38_transcript_cds_lengths.tsv", sep="\t")
cds_data['cds_lengths'] = cds_data['cds_lengths'].apply(eval)
cds_data['cds_exons'] = cds_data['cds_exons'].apply(eval)

clinvar_date = list(Path(str(app_dir) + "/genome_files/clinvar/").glob("*.gz"))[0].name.split("_")[1]
clinvar_date = f"{clinvar_date[0:4]}-{clinvar_date[4:6]}-{clinvar_date[6:8]}"

question_circle_fill = ui.HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle-fill mb-1" viewBox="0 0 16 16"><path d="M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.496 6.033h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286a.237.237 0 0 0 .241.247zm2.325 6.443c.61 0 1.029-.394 1.029-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94 0 .533.425.927 1.01.927z"/></svg>')

# Add page title and sidebar
app_ui = ui.page_navbar(
    ui.head_content(
        ui.HTML("<!-- Google tag (gtag.js) -->\n<script async src='https://www.googletagmanager.com/gtag/js?id=G-MFQGJY5LJY'></script>\n<script>\n  window.dataLayer = window.dataLayer || [];\n  function gtag(){dataLayer.push(arguments);}\n  gtag('js', new Date());\n\n  gtag('config', 'G-MFQGJY5LJY');\n</script>"),
    ),
    ui.nav_panel("epegRNA Designer", 
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
            ui.input_checkbox_group(  
                "checkbox_group",  
                ui.span("Include additional variants in editing windows:", style="text-decoration: underline;"),  
                {  
                    "BLB": "Benign/Likely Benign variants",  
                    "PLP": "Pathogenic/Likely Pathogenic variants",  
                    "GNOMAD": "gnomAD variants with allele count >= 5",
                    "PTC": "PTC variants to induce LoF for assay validation/calibration",  
                },
                width="100%", ),
            ui.span(ui.input_task_button("action_button", "Design epegRNA Libraries"), style="align-self: center;"),
            ui.span(ui.code(f"ClinVar Version: {clinvar_date}"), style="align-self: center;"),
            ui.span(),
            ui.accordion(ui.accordion_panel('Additional Options',
                ui.input_numeric("length_rtt", label="RTT Length", value=40, step=1, min=0, max=60),
                ui.input_numeric("length_pbs", label="PBS Length", value=10, step=1, min=4, max=20),
                ui.input_checkbox("disrupt_pam", "Disrupt PAM/Seed region? (synonymous variants)", value=True, width="100%"),
                ui.input_checkbox("excl_dup_guides", "Exclude spacers duplicated in genome?", value=True, width="100%"),
                ui.input_file("gnomad_csv", ui.tooltip(ui.span("Choose gnomAD csv to upload  ", question_circle_fill),"This file is generated from gnomAD after searching for your gene, filtering as desired, and clicking \"Export variants\"",placement="right"), multiple=False),
                ui.input_numeric("allele_min", label="gnomAD minimum allele count", value=5, step=1, min=0),
                ),id="additional_options", multiple=False, open=False),
            open="desktop",
            width=375,
            padding=20,
            fillable=True,
        ),
        ui.layout_columns(
            ui.value_box(
                "# epegRNA designs", ui.output_ui("total_pegs"),  theme = ui.value_box_theme(bg = "#e6f2fd", fg = "#0B538E")
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
                ui.card_header("epegRNA designs"), ui.output_data_frame("render_design_table"), full_screen=True
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
                    "epegRNA distribution",
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
    ui.nav_control(ui.a('CliPE Home', href='https://home.clipe-mave.org', target="_blank",)),
    ui.nav_control(ui.a('Carvill Lab', href='https://sites.northwestern.edu/carvilllab/', target="_blank",)),
    ui.nav_control(ui.a(ui.img(src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png", alt="Github", height="25px"), href='https://github.com/nbodkin/clipe_pegrna_designer', target="_blank")),

    title=ui.img(src="clipe_logo.png", alt="clipe logo", height="50px"),
    window_title="CliPE epegRNA Designer",
    fillable=True,
)

def server(input, output, session):
    peg_df_glob = reactive.value(pd.DataFrame())
    arch_df_glob = reactive.value(pd.DataFrame())
    fish_df_glob = reactive.value(pd.DataFrame())
    fish_fa_txt_glob = reactive.value("")

    ui.update_selectize("gene", choices=gene_names, selected="TSC2", server=True)

    @reactive.effect
    @reactive.event(input.gene)
    def _():
        if input.gene() == "":
            choices = []
        else:
            choices = gene_data_dict[input.gene()]
            ui.update_select("transcript", choices=choices)
        
    
    @render.data_frame
    def render_design_table():
        return async_run_guide_design.result()

    @reactive.effect
    @reactive.event(input.action_button, ignore_init=True)
    def build_peg_df():
        async_run_guide_design(input.gene(), input.transcript(), input.length_pbs(), input.length_rtt(), input.num_designs(), input.design_strategy(), input.checkbox_group(), input.allele_min(), input.gnomad_csv(), input.disrupt_pam(), input.excl_dup_guides())
        
    
    @ui.bind_task_button(button_id="action_button")
    @reactive.extended_task
    async def async_run_guide_design(gene_name, transcript, length_pbs, length_rtt, num_designs, design_strategy, checkbox_group, allele_min, gnomad_csv, disrupt_pam, excl_dup_guides):
        # remove download button if it exists
        ui.remove_ui("#download_button")
        ui.remove_ui("#download_checkbox")
        peg_df_glob.set(pd.DataFrame())
        
        # validate input files exist
        if not gnomad_csv:
            gnomad_file = None
            if "GNOMAD" in checkbox_group:
                ui.notification_show("Please upload a gnomAD file to include gnomAD variants in the design.", duration=5, type="warning")
        else:
            gnomad_file = gnomad_csv[0]["datapath"]
        
        with ui.Progress(min=0, max=6) as p:
            p.set(0, message="Running")
            expt = clipe_expt(gene_name, transcript.split(" ")[0], gnomad_file, length_pbs, length_rtt, num_designs, disrupt_pam, excl_dup_guides, design_strategy, checkbox_group, allele_min, prog_bar=p)
            peg_df, arch_df, windows = expt.run_guide_design()
            p.set(5, detail="Finishing up")
            peg_df_glob.set(peg_df)
            arch_df_glob.set(arch_df)

            output = build_files_for_jellyfish(peg_df, screening_df=arch_df)
            fish_df_glob.set(output[0])
            fish_fa_txt_glob.set(output[1])
            
            #free up memory
            del expt
            gc.collect()
            p.set(6)
        
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
            "peg_tables": "epegRNA designs and ordering (.csv, .xlsx)", 
            "arch_tables": "archetypal epegRNA files (.csv, .txt)",
            "RTTs": "RTT data for downstream analyses (.csv, .fa)",
            "nicking_idt": "IDT nicking guide ordering files (.txt)"
            },
            width="100%", selected=["peg_tables", "arch_tables", "RTTs", "nicking_idt"]),
            selector="#download_area",
            where="afterEnd",
        )
        
        return peg_df

    @render.download(filename=lambda: f"{date.today()}_{input.gene()}_clipe_designs.zip")
    def download_button():
        peg_df = peg_df_glob.get()
        arch_df = arch_df_glob.get()
        fish_df = fish_df_glob.get()
        fish_fa_txt = fish_fa_txt_glob.get()
        files_to_download = input.download_checkbox()
        if len(files_to_download) > 0:
            file_prefix = f"{input.gene()}_"
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                if "peg_tables" in files_to_download:
                    peg_df.to_csv(temp_path / f"{file_prefix}full_epegRNA_designs.csv", index=False)
                    idt_df = peg_df[['editing_window', "full_peg"]]
                    idt_df['editing_window'] = idt_df['editing_window'].apply(lambda x: f"{file_prefix}window_{x}")
                    idt_df.columns = ["Pool name", "Sequence"]
                    idt_df.to_excel(temp_path / f"{file_prefix}full_IDT_opool_order_data.xlsx", index=False)
                if "arch_tables" in files_to_download:
                    arch_df.to_csv(temp_path / f"{file_prefix}archetypal_epegRNA_designs.csv", index=False)
                    idt_df = prep_screening_order_df(file_prefix, arch_df)
                    idt_df[['name', 'oligo', 'scale', 'purification']].to_csv(temp_path / f"{file_prefix}archetypal_IDT_oligo_order_data.txt", sep="\t", index=False, header=False)
                if "RTTs" in files_to_download:
                    with open(temp_path / f"{file_prefix}RTT_fasta.fa", "w") as f:
                        f.write(fish_fa_txt)
                    fish_df.to_csv(temp_path / f"{file_prefix}RTT_table.csv", index=False)                    
                if "nicking_idt" in files_to_download:
                    nick_df = prep_nicking_order_df(file_prefix, peg_df)
                    nick_df[['name', 'oligo', 'scale', 'purification']].to_csv(temp_path / f"{file_prefix}nicking_IDT_oligo_order_data.txt", sep="\t", index=False, header=False)
                    
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
            return peg_df[peg_df['Germline classification'].isin(["Uncertain_significance", "Conflicting_classifications_of_pathogenicity"])].shape[0]
    
    @render.ui
    def total_plp():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])].shape[0]
    
    @render.ui
    def total_blb():
        peg_df = peg_df_glob.get()
        if peg_df.shape[0] > 0:
            return peg_df[peg_df['Germline classification'].isin(["Benign", "Likely_benign", "Benign/Likely_benign"])].shape[0]

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