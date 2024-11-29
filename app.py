from datetime import date
from clipe_guide_design_methods import *
from shiny import App, render, ui, reactive
from pathlib import Path


# A card component wrapper.
def ui_card(title, *args):
    return (
        ui.div(
            {"class": "mx-auto card mb-4 col-sm-10 col-md-8 col-lg-6 "},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body"}, *args),
        ),
    )

app_ui = ui.page_fluid(
    ui.h2(
        "CliPE: Curated Loci Prime Editing: pegRNA Design Tool",
        class_=" mx-auto text-lg-center text-left py-5",
    ),
    ui.div(class_="pb-4"),
    ui_card(
        "CliPE Parameters",
        ui.input_radio_buttons(
         "design_strategy",
         "Choose an experiment design strategy:",
         {
             "vus": "Identify regions with highest denisty of clinvar missense VUS",
             "pos_control": "Identify regions with highest denisty of Clinvar missense BLB and PLP",
         },
         width="100%",
        ),
        ui.input_text("GeneName", label="Gene Name", value="TSC2"),
        ui.input_text("transcript", label="Transcript (Refseq)", value="NM_000548.5"),
        ui.input_numeric("num_designs", label="Number of epegRNA libraries to design", value=12, step=1),
        ui.input_numeric("length_rtt", label="Length of RTT", value=40, step=1),
        ui.input_numeric("length_pbs", label="Length of PBS", value=10, step=1),
        ui.input_radio_buttons("disrupt_pam", "Disrupt PAM/Seed region with synonymous variants?", {"yes": "Yes", "no": "No"}, width="100%"),
        ui.input_file("clinvar_CSV", "Choose missense Clinvar csv to upload:", multiple=False),
        ui.input_file("gnomAD_CSV", "Choose gnomAD csv to upload:", multiple=False),
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
        ui.input_action_button("action_button", "Design epegRNA Libraries")
        ),

    ui_card(
        "pegRNA Libraries",
        ui.output_data_frame("build_peg_df"),
    ),

)


def server(input, output, session):
    @render.data_frame
    @reactive.event(input.action_button)
    def build_peg_df():
        # remove download button if it exists
        ui.remove_ui("#download")

        # validate input files exist
        if not input.clinvar_CSV():
            raise ValueError("Clinvar file must be provided")
        if not input.gnomAD_CSV():
            gnomad_file = None
        else:
            gnomad_file = input.gnomAD_CSV()[0]["datapath"]
        
        if input.disrupt_pam() == "yes":
            disrupt_pam = True
        else:
            disrupt_pam = False
    
        expt = clipe_expt(input.transcript(), input.clinvar_CSV()[0]["datapath"], gnomad_file, input.length_pbs(), input.length_rtt(), input.num_designs(), disrupt_pam, input.design_strategy(), input.checkbox_group(), input.allele_min())
        global peg_df
        peg_df = expt.run_guide_design()

        ui.insert_ui(
            ui.download_button("download", "Download CSV"),
            selector="#build_peg_df",
            where="afterEnd",
        )

        return peg_df

    @render.download()
    def download():
        path = Path(__file__).parent / f"{input.GeneName()}_clipe_peg_df_{date.today()}.csv"
        print(path)
        peg_df.to_csv(path, index=False)
        return str(path)


app = App(app_ui, server)
