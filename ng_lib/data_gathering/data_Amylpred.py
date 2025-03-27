import time
import requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
import os
import re
# Define the CDC19 protein sequence in FASTA format
protein_fasta_full = """>sp|P00549|KPYK1_YEAST Pyruvate kinase 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=CDC19 PE=1 SV=2
MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYH
KSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKA
CDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVN
LPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVK
IENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQM
LESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIA
YLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPI
ILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYV
SIQGFKAGAGHSNTLQVSTV"""

AMYLPRED_DEBUG = True


def fetch_amylpred_results(fasta_sequence, username, password):
    url = "http://thalis.biol.uoa.gr/AMYLPRED2/input.php"
    driver = webdriver.Chrome()
    driver.get(url)

    textarea = driver.find_element(By.NAME, "email")
    textarea.send_keys(username)
    
    textarea = driver.find_element(By.NAME, "password")
    textarea.send_keys(password)
    
    # Find the textarea and input the sequence
    login_button = driver.find_element(By.XPATH, '//input[@type="submit" and @value="Login"]')
    login_button.click()

    textarea = driver.find_element(By.NAME, "seq_data")
    textarea.send_keys(fasta_sequence)
    
    submit_button = driver.find_element(By.XPATH, '//input[@class="normalc" and @value="Submit Query" and @type="submit"]')
    submit_button.click()

    driver.implicitly_wait(120)    
    link = driver.find_element(By.XPATH, '//a[contains(@href, "/AMYLPRED2/tmp/") and contains(@href, ".txt")]')
    link.click()

    # Extract the URL from the href attribute
    results_url = link.get_attribute('href')

    driver.quit()
    
    return results_url

def download_results_file(url):
    # TODO: just return response.content
    response = requests.get(url)
    filename = "Amylpred_CDC19_results.txt"
    with open(filename, 'wb') as file:
        file.write(response.content)
    return filename

def parse_results_file(filename):
    # Read the content of the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    hit_section = False
    consensus_dict = {}

    # Process each line
    for line in lines:
        # Check if we've reached the HITS section
        if "HITS" in line:
            hit_section = True
            continue

        if hit_section:
            # Match lines that start with CONSENSUS
            line = line.replace(">--->", "     ")
            match = re.match(r'\s*(CONSENSUS\d+):\s*(.*)', line)
            if match:
                consensus_type = match.group(1)
                ranges = match.group(2).strip()
                # If there are ranges, split them into a list
                if ranges:
                    range_list = ranges.split(', ')
                    consensus_dict[consensus_type] = range_list
                else:
                    consensus_dict[consensus_type] = []
        
    return consensus_dict

def get_consensus_vec(consensus_dict, protein_length):
    cons = "CONSENSUS"
    consensus_dict_out = dict()
    consensus_vec = np.zeros(protein_length)
    for i in range(2,10):
        for r in consensus_dict[cons+str(i)]:
            r_doubles = r.split('-')
            range_start = int(r_doubles[0])
            range_end = int(r_doubles[1])
            consensus_vec[range_start:range_end] +=1
    return consensus_vec

def get_Amylpred_data(protein_fasta, plot_it = False):

    if AMYLPRED_DEBUG:
        print("data_Amylpred: DEBUG MODE: using CDC19 data")
        # Get the directory of the current script
        current_script_dir = os.path.dirname(__file__)

        # Construct the path to the file
        results_filename = os.path.join(current_script_dir, "Amylpred_CDC19_results.txt") 
    else:
            
        # Fetch the results URL
        results_url = fetch_amylpred_results(protein_fasta)

        # Download the results file

        results_filename = download_results_file(results_url)

    # Parse the results file
    consensus_dict = parse_results_file(results_filename)
    consensus_vec = get_consensus_vec(consensus_dict, len(protein_fasta))

    # Create a list to hold the scores
    sequence_only = protein_fasta
    protein_length = len(sequence_only)
    scores = consensus_vec.reshape(1, -1)


    plot_it = True
    if plot_it:
        # Plotting the results as a heatmap
        plt.figure(figsize=(10, 2))
        plt.imshow(scores, cmap='viridis', aspect='auto')
        plt.colorbar(label='Score')
        plt.yticks([])
        plt.xticks(ticks=np.arange(protein_length), fontsize=8)
        plt.xlabel('Amino Acid Position')
        plt.title('Protein Sequence Analysis - CDC19')
        plt.show()
    return scores


def gather_cdc19_data():
    cdc19Fasta = "MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYH\
    KSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKA\
    CDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVN\
    LPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVK\
    IENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQM\
    LESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIA\
    YLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPI\
    ILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYV\
    SIQGFKAGAGHSNTLQVSTV"
    get_zipperDB_data(cdc19Fasta, plot_it = True)
    


if __name__ == "__main__":
    get_Amylpred_data(protein_fasta_full)