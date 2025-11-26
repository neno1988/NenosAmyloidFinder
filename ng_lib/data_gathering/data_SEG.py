import time
import re
import requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import os
import filecache

SEG_DEBUG = False

@filecache.filecache(7*24*60*60)
def fetch_seg_results(fasta_sequence):
    if SEG_DEBUG:
        print("data_SEG: DEBUG MODE: using CDC19 data")
        # Get the directory of the current script
        current_script_dir = os.path.dirname(__file__)

        # Construct the path to the file
        file_path = os.path.join(current_script_dir, 'CDC19_SEG_result.txt')
        with open(file_path, "r") as file:
            txt = file.read() 
        return txt
    
    options = Options()
    #options.add_argument("--headless")  # Run Chrome in headless mode
    options.add_argument("--disable-gpu")  # Optional: avoid GPU-related issues on Windows
    options.add_argument("--window-size=1920,1080")  # Optional: simulate full screen for page rendering
    driver = webdriver.Chrome(options = options)

    url = "https://mendel.imp.ac.at/METHODS/seg.server.html"
    driver.get(url)

    # Find the textarea and input the sequence
    textarea = driver.find_element(By.NAME, "Sequence")
    textarea.send_keys(fasta_sequence.seq)
    
    # Find the submit button and click it
    submit_button = driver.find_element(By.XPATH, "//input[@type='SUBMIT']")
    submit_button.click()
    
    # Wait for the results to load
    time.sleep(8)  # Adjust the time based on the actual processing time

    # Find the results element and get the text
    results_element = driver.find_element(By.XPATH, "//pre")
    results_text = results_element.text

    driver.quit()
    
    return results_text

def parse_seg_results(results_text):
    # Parse the SEG results from the text output
    segments_12 = []
    segments_25 = []
    segments_45 = []
    
    segments_dict = {'SEG 12': segments_12, 'SEG 25': segments_25, 'SEG 45': segments_45}
    current_segment = None
    
    pattern = re.compile(r"(\d+)-(\d+)")
    pattern_low_complexity = re.compile(r"([A-Za-z]+)\s+(\d+)-(\d+)\s*")
    pattern_not_low_complexity = re.compile(r"\s*(\d+)-(\d+)\s+([A-Za-z]+)")
    
    for line in results_text.split('\n'):
        if 'SEG 12' in line:
            current_segment = 'SEG 12'
        elif 'SEG 25' in line:
            current_segment = 'SEG 25'
        elif 'SEG 45' in line:
            current_segment = 'SEG 45'
        elif current_segment and pattern.search(line):
            match = pattern.search(line)
            match_LC = pattern_low_complexity.search(line) 
            match_not_LC = pattern_not_low_complexity.search(line)
            assert(match_LC != match_not_LC)
            if match:
                start, end = match.groups()
                start, end = int(start), int(end)
                segments_dict[current_segment].append((start, end, bool(match_LC)))
                
    return segments_12, segments_25, segments_45

def get_SEG_data(protein_fasta, plot_it = False, debug = False, seg=12):
    # Fetch and parse the results
    seg_results_text = fetch_seg_results(protein_fasta)
    seg_results_12, seg_results_25, seg_results_45 = parse_seg_results(seg_results_text)

    # Create a list to hold the features
    sequence_only = protein_fasta.seq# "".join(protein_fasta.split('\n')[2:])
    protein_length = len(sequence_only)
    features = np.zeros(protein_length)

    # Use the SEG 12 results for the heatmap
    if seg == 0:
        seg_results = seg_results_12 #TODO: hack to be solved
    elif seg==12:
        seg_results = seg_results_12
    elif seg==25:
        seg_results = seg_results_25
    elif seg==45:
        seg_results = seg_results_45
    else:
        raise ValueError("Invalid SEG number")
    
    for start, end_seq, isLC in seg_results:
        features[start-1:end_seq] = float(isLC) 

    # Plotting the results as a heatmap
    if plot_it:
        plt.figure(figsize=(10, 2))
        plt.imshow(features.reshape(1, -1), cmap='viridis', aspect='auto')
        plt.colorbar(label='Feature')
        plt.yticks([])
        plt.xticks(ticks=np.arange(protein_length), labels=list(sequence_only), fontsize=8)
        plt.xlabel('Amino Acid Position')
        plt.title('Protein Sequence Analysis - CDC19 (SEG 12)')
        plt.show()
    return features


# TODO: use shared class, not this custom made one.

class FastaSeq:
    seq:str
    description:str

    def __init__(self, seq="", description=""):
        self.seq = seq
        self.description = description


def gather_cdc19_data():

    # Define the CDC19 protein sequence in FASTA format
    cdc19Fasta = FastaSeq()
    cdc19Fasta.seq = "MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYHKSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKACDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVNLPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVKIENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQMLESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIAYLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPIILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYVSIQGFKAGAGHSNTLQVSTV"
    cdc19Fasta.description = "sp Some description (no > at the beginning)"
    
    get_SEG_data(cdc19Fasta, plot_it = True)


if __name__ == "__main__":
    gather_cdc19_data()