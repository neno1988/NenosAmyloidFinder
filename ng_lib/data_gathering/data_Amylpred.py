import time
import requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import os
import re
import filecache
import json
from .interface import DataGatheringTool

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

CDC19_AMYLPRED2_RESULT_FILE =r"\ng_lib\data_gathering\CDC19_Amylpred2_results.txt"

AMYLPRED_DEBUG = False
DEBUG_SECRET_FILE_PATH = r"C:\Users\Neno\Documents\Git\NenoGeaProject\NenosAmyloidFinder\neno_config.json"



class AmylpredDataGatheringTool(DataGatheringTool):
    def get_data_from_sequence(self, sequence) -> np.ndarray:
        return self.get_amylpred_data(sequence)


    @filecache.filecache(7*24*60*60)
    def fetch_amylpred_results(self, fasta_sequence):
        options = Options()
        #options.add_argument("--headless")  # Run Chrome in headless mode
        options.add_argument("--disable-gpu")  # Optional: avoid GPU-related issues on Windows
        options.add_argument("--window-size=1920,1080")  # Optional: simulate full screen for page rendering
        driver = webdriver.Chrome(options = options)

        url = "http://thalis.biol.uoa.gr/AMYLPRED2/input.php"
        driver.get(url)

        # Get credentials from environment variables
        username = os.environ.get('AMYLPRED_USERNAME') # TODO: get from GUI is not available in environment
        password = os.environ.get('AMYLPRED_PASSWORD') # TODO: get from GUI is not available in environment

        if not username or not password:
            if os.path.exists(DEBUG_SECRET_FILE_PATH): # TODO: create and use configuration handler class or functions
                with open(DEBUG_SECRET_FILE_PATH, 'r') as f:
                    config = json.load(f)
            username = config.get("amylpred_username")
            password = config.get("amylpred_password")
            
        if not username or not password:
            raise ValueError("Amylpred credentials not found. Please provide username and password in the GUI.")

        textarea = driver.find_element(By.NAME, "email")
        textarea.send_keys(username)
        
        textarea = driver.find_element(By.NAME, "password")
        textarea.send_keys(password)
        
        # is it a batch sequence?
        if fasta_sequence.count('>')>1:
            print("Batch mode! Not ready!")
            raise(ValueError("Batch mode not available for Amylpred"))

        else:
            # Find the textarea and input the sequence
            login_button = driver.find_element(By.XPATH, '//input[@type="submit" and @value="Login"]')
            login_button.click()

            textarea = driver.find_element(By.NAME, "seq_data")
            textarea.send_keys(fasta_sequence)
            
            #uncheck AMYLPATTERN and 
            # Find the checkbox by name, value, or another attribute
            checkbox_aggrescan = driver.find_element(By.XPATH, "//input[@type='checkbox' and @name='method' and @value='AGGRESCAN']")
            checkbox_amylmuts = driver.find_element(By.XPATH, "//input[@type='checkbox' and @name='method' and @value='AMYLMUTS']")
            # Check if it is selected
            if checkbox_aggrescan.is_selected():
                checkbox_aggrescan.click()  # uncheck the box
            if checkbox_amylmuts.is_selected():
                checkbox_amylmuts.click()  # uncheck the box

            submit_button = driver.find_element(By.XPATH, '//input[@class="normalc" and @value="Submit Query" and @type="submit"]')
            submit_button.click()

            driver.implicitly_wait(120)    
            link = driver.find_element(By.XPATH, '//a[contains(@href, "/AMYLPRED2/tmp/") and contains(@href, ".txt")]')
            link.click()

            # Extract the URL from the href attribute
            results_url = link.get_attribute('href')

            driver.quit()
            
        return results_url


    def download_results_file(self, url):
        # TODO: just return response.content
        response = requests.get(url)
        filename = "Amylpred_CDC19_results.txt"
        with open(filename, 'wb') as file:
            file.write(response.content)
        return filename

    def parse_results_file(self, filename):
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

    def get_consensus_vec(self, consensus_dict, protein_length):
        cons = "CONSENSUS"
        consensus_dict_out = dict()
        consensus_vec = np.zeros(protein_length)
        for i in range(2,8):
            for r in consensus_dict[cons+str(i)]:
                r_doubles = r.split('-')
                range_start = int(r_doubles[0])
                range_end = int(r_doubles[1])
                consensus_vec[range_start:range_end] +=1
        return consensus_vec

    def get_amylpred_data(self, protein_fasta, plot_it = False):

        if AMYLPRED_DEBUG:
            print("data_Amylpred: DEBUG MODE: using CDC19 data")
            # Get the directory of the current script
            current_script_dir = os.path.dirname(__file__)

            # Construct the path to the file
            results_filename = os.path.join(current_script_dir, "Amylpred_CDC19_results.txt") 
        else:
            results_url = self.fetch_amylpred_results(protein_fasta)
            results_filename = self.download_results_file(results_url)

        # Parse the results file
        consensus_dict = self.parse_results_file(results_filename)
        consensus_vec = self.get_consensus_vec(consensus_dict, len(protein_fasta))

        # Create a list to hold the scores
        sequence_only = protein_fasta
        protein_length = len(sequence_only)
        scores = consensus_vec.reshape(1, -1)


        plot_it = False
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

def test_single_fasta():
    prepare_environment()
    tool = AmylpredDataGatheringTool()
    res_url = tool.fetch_amylpred_results(CDC19_FASTA_FULL)
    assert(res_url[0:35]=="http://thalis.biol.uoa.gr/AMYLPRED2")
    print(res_url)


def test_multi_fasta():
    tool = AmylpredDataGatheringTool()
    prepare_environment()
    tool.fetch_amylpred_results(CDC19_FASTA_FULL + "\n\n" + CDC19_FASTA_FULL)


def test_consensus_vec():
    tool = AmylpredDataGatheringTool()
    test_result_file_path = os.getcwd() + CDC19_AMYLPRED2_RESULT_FILE
    with open(test_result_file_path, 'r') as file:
        lines = file.readlines()
    consensus_dict = tool.parse_results_file(test_result_file_path)
    consensus_vec = tool.get_consensus_vec(consensus_dict, len(CDC19_FASTA_FULL))
    print(consensus_vec)



def load_user_and_pass():
    import json
    SECRET_FILE = "secret.json"
    if os.path.exists(SECRET_FILE):
        with open(SECRET_FILE, 'r') as f:
            return json.load(f)
    return {}

def prepare_environment():
        user_pass_dict = load_user_and_pass()
            # Set Amylpred credentials in the environment
        os.environ['AMYLPRED_USERNAME'] = user_pass_dict["amylpred_username"]
        os.environ['AMYLPRED_PASSWORD'] = user_pass_dict["amylpred_password"]

if __name__ == "__main__":
    # get_amylpred_data(protein_fasta_full)
    # test_single_fasta()
    test_consensus_vec()
