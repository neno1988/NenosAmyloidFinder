import sys
from pathlib import Path

# add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import plotly as plt
import ng_lib.data_gathering as dg
from ng_lib.utils import FastaSeq
from ng_lib.plot import heatmaps_binary_non_binary
from ng_lib.visualization.plotly_heatmaps import create_interactive_heatmaps, make_annotations, HeatmapElement
import os
import json
import re
import numpy as np


# Define the CDC19 protein sequence in FASTA format
CDC19_FASTA_FULL = """>sp|P00549|KPYK1_YEAST Pyruvate kinase 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=CDC19 PE=1 SV=2
MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYH
KSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKA
CDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVN
LPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVK
IENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQM
LESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIA
YLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPI
ILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYV
SIQGFKAGAGHSNTLQVSTV"""


def analyse_protein(output_folder, seq, description, name, threshold, SEG, xticks, DEBUG = False):
    #TODO: this function might belong in the upper sub-package
    fasta_protein = FastaSeq() # TODO: use official fasta struct here
    fasta_protein.description = description
    fasta_protein.seq = seq
    fasta_protein.name = name

    # Get SEG and ZipperDB data
    seg_tool_parameters = dg.SEGDataGatheringToolParameters(seg=SEG)
    seg_tool = dg.SEGDataGatheringTool()
    if DEBUG:
        seg_data = seg_tool.get_debug_data(len(fasta_protein.seq))
    else:
        seg_data = seg_tool.get_data_from_sequence(fasta_protein, parameters=seg_tool_parameters)
    zipperDB_parameters = dg.ZDBDataGatheringToolParameters(threshold=threshold)
    zipperDB_tool = dg.ZDBDataGatheringTool()

    if DEBUG:
        zipperDB_data = zipperDB_tool.get_debug_data(len(fasta_protein.seq))
    else:
        zipperDB_data = zipperDB_tool.get_data_from_sequence(fasta_protein, parameters=zipperDB_parameters)

    amylpred_tool = dg.AmylpredDataGatheringTool()
    # Get Amylpred data if credentials are provided
    NUMBER_OF_AMYLPRED_TOOLS = 6
    if DEBUG:
        amylpred_data = amylpred_tool.get_debug_data(len(fasta_protein.seq))*6
    else:
        try:
            amylpred_data = amylpred_tool.get_data_from_sequence(fasta_protein.seq)
        except Exception as e:
            # TODO: reinstate warning
            # messagebox.showwarning("Warning", f"Failed to get Amylpred data: {str(e)}")
            print("Reinstate warning about Amylpred failure here")
            amylpred_data = None

    aggrescan_tool = dg.AggrescanDataGatheringTool()
    if DEBUG:
        aggrescan_data = aggrescan_tool.get_debug_data(len(fasta_protein.seq))
    else:
        aggrescan_data = aggrescan_tool.get_data_from_sequence(fasta_protein.seq)
    
    # complete amylopred data with the typically missing aggrescan data
    amylpred_data = aggrescan_data + (amylpred_data if amylpred_data is not None else np.zeros_like(aggrescan_data))
    NUMBER_OF_AMYLPRED_TOOLS += 1
    
    # Create interactive plot
    IMAGE_HEIGHT = 200
    IMAGE_WIDTH = 1800
    heatmaps: list[HeatmapElement] = []
    heatmaps.append(HeatmapElement(seg_data, number_of_rows=2, min_value=0, max_value=1, colors = ["White", "orange"], legend = "LCR"))
    heatmaps.append(HeatmapElement(amylpred_data, number_of_rows=7, min_value=0, max_value=NUMBER_OF_AMYLPRED_TOOLS, colors = ["White", "green"], legend = "Amylpred"))
    heatmaps.append(HeatmapElement(zipperDB_data, number_of_rows=2, min_value=0, max_value=1, colors = ["White", "darkgray"], legend = "ZipperDB"))
 

    fig = create_interactive_heatmaps(
        heatmaps=heatmaps,
        name=name,
        xticks=xticks,
        image_height=IMAGE_HEIGHT, 
        image_width=IMAGE_WIDTH 
    )
    make_annotations(fig, name, "left", SEG, threshold)
    #fig.show()

    # Save the figure as HTML
    name = re.sub(r'[^\w_. -]', '_', name)
    name_and_parameters = f"{name}_SEG{SEG}_ZDB{threshold}"
    #output_folder = os.path.dirname(output_folder)
    html_file = os.path.join(output_folder, f"{name_and_parameters}_analysis.html")
    html_file = os.path.normpath(html_file)
    png_file = os.path.join(output_folder, f"{name_and_parameters}_image.png")
    json_file = os.path.join(output_folder, f"{name_and_parameters}_data.json")

    fig.write_html(html_file)
    fig.data= [fig.data[0]]
    fig.update_layout(title_text=None)
    fig.write_image(png_file, engine='kaleido', width = IMAGE_WIDTH, height = IMAGE_HEIGHT)


    # Generate output data file
    data_file = dict()
    data_file["fasta_seq"] = fasta_protein.seq
    data_file["fasta_description"] = fasta_protein.description
    data_file["fasta_name"] = fasta_protein.name
    data_file["ZipperDB"] = zipperDB_data.tolist()
    data_file["Amylpred"] = amylpred_data.tolist()
    data_file["SEG"] = seg_data.tolist()

    # Save the JSON
    with open(json_file, 'w') as f:
        json.dump(data_file, f, indent = 4)
    
    return fig, output_folder


def test_analysis():
    with open(os.getcwd() + f'\\test_files\\test_data_CDC19.json', 'r') as f:
        test_data = json.load(f)
    threshold = -23
    SEG = 25
    xticks = 0
    # TODO: separation of concerns could be better here
    sequence = test_data["fasta_seq"]
    output_folder = os.getcwd() + f'\\test_files\\'
    fig, output_folder = analyse_protein(output_folder, sequence, test_data["fasta_description"], test_data["fasta_name"], threshold, SEG, xticks, DEBUG = True)
    fig.show()

# TEST DATA GENERATION FUNCTIONS
def generate_test_data():
    import numpy as np
    fasta_cdc19 = get_CDC19_fasta()
    SEG_data = np.random.random(len(fasta_cdc19.seq))
    zipperDB_data = np.random.random(len(fasta_cdc19.seq))
    amylpred_data = np.random.random(len(fasta_cdc19.seq))*6

    test_data = {
        "fasta_seq": fasta_cdc19.seq,
        "fasta_description": fasta_cdc19.description,
        "fasta_name": fasta_cdc19.name,
        "SEG_data": SEG_data.tolist(),
        "zipperDB_data": zipperDB_data.tolist(),
        "amylpred_data": amylpred_data.tolist()
    }

    with open(os.getcwd() + f'\\test_files\\test_data_CDC19.json', 'w') as f:
        json.dump(test_data, f, indent=4)
       

def get_CDC19_fasta():
    fasta = FastaSeq()
    fasta.description = """>sp|P00549|KPYK1_YEAST Pyruvate kinase 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=CDC19 PE=1 SV=2"""
    fasta.seq = """MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYH\
KSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKA\
CDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVN\
LPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVK\
IENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQM\
LESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIA\
YLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPI\
ILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYV\
SIQGFKAGAGHSNTLQVSTV"""
    fasta.name = "CDC19"
    return fasta



if __name__ == "__main__":
    # generate_test_data() # generates test data in test_files based on CDC19 fasta.
    test_analysis()
    print("DONE")