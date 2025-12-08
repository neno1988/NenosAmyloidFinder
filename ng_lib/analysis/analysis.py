import sys
from pathlib import Path

# add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import plotly as plt
from ng_lib.data_gathering import get_SEG_data, get_zipperDB_data, get_amylpred_data
from ng_lib.utils import FastaSeq
from ng_lib.plot import heatmaps_binary_non_binary
from ng_lib.visualization.plotly_heatmaps import create_interactive_heatmaps
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


def analyse_protein(output_folder, seq, description, name, threshold, SEG, xticks, 
                    debug_SEG_data = None, debug_zipperDB_data = None, debug_amylpred_data = None):
    #TODO: this function might belong in the upper sub-package
    fasta_protein = FastaSeq() # TODO: use official fasta struct here
    fasta_protein.description = description
    fasta_protein.seq = seq
    fasta_protein.name = name

    # Get SEG and ZipperDB data
    if debug_SEG_data is not None:
        seg_data = debug_SEG_data
    else:
        seg_data = get_SEG_data(fasta_protein, seg=SEG)

    if debug_zipperDB_data is not None:
        zipperDB_data = debug_zipperDB_data
    else:
        zipperDB_data = get_zipperDB_data(fasta_protein)

    # Get Amylpred data if credentials are provided
    if debug_amylpred_data is not None:
        amylpred_data = debug_amylpred_data
    else:
        try:
            amylpred_data = get_amylpred_data(fasta_protein.seq)
        except Exception as e:
            # TODO: reinstate warning
            # messagebox.showwarning("Warning", f"Failed to get Amylpred data: {str(e)}")
            print("Reinstate warning about Amylpred failure here")
            amylpred_data = None

    do_matplotlib_plots = False
    if do_matplotlib_plots:
        if SEG==0:
            # TODO: only return fig  in nice_heatmap_plot?
            white_cmap = ["white", "white"]
            fig = heatmaps_binary_non_binary(seg_data, zipperDB_data[0], threshold=threshold, name=name, force_cmap=white_cmap, seg_bar_height=0, xticks=xticks)
            # Save the figure as a JPEG file
            file_name = fasta_protein.name+'_treshold_'+str(threshold)+'.jpg'
            fig.savefig(os.path.join(output_folder, file_name), format='jpeg', dpi=300)
            fig.show() 

        else:
            fig = heatmaps_binary_non_binary(seg_data, zipperDB_data[0], threshold=threshold, name=name, xticks=xticks)
            # Save the figure as a JPEG file
            file_name = fasta_protein.name+'_treshold_'+str(threshold)+'.jpg'
            fig.savefig(os.path.join(output_folder, file_name), format='jpeg', dpi=300)
            fig.show() 
    
    # Create interactive plot
    IMAGE_HEIGHT = 200
    IMAGE_WIDTH = 1800
    fig = create_interactive_heatmaps(
        lcr_data=seg_data.reshape(1, -1),
        seg_treshold = SEG,
        zipperdb_data=zipperDB_data,
        amylpred_data=amylpred_data if amylpred_data is not None else np.zeros_like(zipperDB_data),
        zipperdb_threshold=threshold,
        name=name,
        xticks=xticks,
        image_height=IMAGE_HEIGHT, 
        image_width=IMAGE_WIDTH 
    )
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
    
    # If Amylpred data is available, show it
    plot_amylpred_only = False
    if plot_amylpred_only and amylpred_data is not None:
        plt.figure(figsize=(10, 2))
        plt.imshow(amylpred_data, cmap='viridis', aspect='auto')
        plt.colorbar(label='Amylpred Score')
        plt.yticks([])
        plt.xticks(ticks=np.arange(len(seq)), fontsize=8)
        plt.xlabel('Amino Acid Position')
        plt.title(f'Amylpred Analysis - {name}')
        plt.show()


def test_analysis():
    with open(os.getcwd() + f'\\test_files\\test_data_CDC19.json', 'r') as f:
        test_data = json.load(f)
    threshold = -23
    SEG = 25
    xticks = 0
    # TODO: separation of concerns could be better here
    sequence = test_data["fasta_seq"]
    output_folder = os.getcwd() + f'\\test_files\\'
    analyse_protein(output_folder, sequence, test_data["fasta_description"], test_data["fasta_name"], threshold, SEG, xticks, 
                    debug_SEG_data = np.array(test_data["SEG_data"]), 
                    debug_zipperDB_data = np.array(test_data["zipperDB_data"]), 
                    debug_amylpred_data = np.array(test_data["amylpred_data"]))


# TEST DATA GENERATION FUNCTIONS
def generate_test_data():
    fasta_cdc19 = get_CDC19_fasta()
    SEG_data = get_SEG_data(fasta_cdc19, seg=25)
    zipperDB_data = get_zipperDB_data(fasta_cdc19)
    amylpred_data = get_amylpred_data(fasta_cdc19.seq)

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