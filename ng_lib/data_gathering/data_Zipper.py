import time
import requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
import os
import filecache

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

ZIPPERDB_DEBUG = False

def fetch_single_zipperdb_result(fasta_sequence):
    url = "https://zipperdb.mbi.ucla.edu"
    driver = webdriver.Chrome()
    driver.get(url)

    # Find the textarea and input the sequence
    textarea = driver.find_element(By.ID, "tfasta")
    textarea.send_keys(fasta_sequence)

    # Find the other input and submit (assuming it's needed)
    pname = driver.find_element(By.ID, "pname")
    pname.send_keys("a")

    # Assuming there's a submit button and its ID is 'submit_btn'
    submit_btn = driver.find_element(By.ID, "submit_btn")
    submit_btn.click()

    # Wait for the results to load (you might need to add explicit waits here)
    driver.implicitly_wait(10)

    # Find the table containing the results
    table = driver.find_element(By.CSS_SELECTOR, "div.tableBox table.grid")

    # Extract table headers
    headers = [header.text for header in table.find_elements(By.TAG_NAME, "th")]

    # Extract table rows
    rows = table.find_elements(By.TAG_NAME, "tr")
    data = []
    for row in rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        if cells:
            data.append([cell.text for cell in cells])

    # Close the driver
    driver.quit()

    # Create a pandas DataFrame
    df = pd.DataFrame(data, columns=headers)
    return df["Score"].reshape(1,-1)
    
    


def fetch_zipperdb_results(fasta_sequence):
    url = "https://zipperdb.mbi.ucla.edu/batch"
    driver = webdriver.Chrome()
    driver.get(url)

    # Find the textarea and input the sequence
    textarea = driver.find_element(By.ID, "tfasta")
    textarea.send_keys(">")
    textarea.send_keys(fasta_sequence.description)
    textarea.send_keys("\n")
    textarea.send_keys(fasta_sequence.seq)
    
    # Find the submit button and click it
    submit_button = driver.find_element(By.ID, "dobatch")
    submit_button.click()
    
    # Wait for the results link to appear
    time.sleep(6)  # Adjust the time based on the actual processing time

    # Find the results link
    results_link = driver.find_element(By.XPATH, "//a[contains(text(), 'Collect results here')]")
    results_url = results_link.get_attribute('href')

    driver.quit()
    
    return results_url

def download_results_file(url):
    # TODO: just return response.content
    response = requests.get(url)
    filename = "zipperdb_results.csv"
    with open(filename, 'wb') as file:
        file.write(response.content)
    return filename

def parse_results_file(filename):
    df = pd.read_csv(filename, comment='#')
    return df

@filecache.filecache(24*60*60)
def get_zipperDB_data(protein_fasta, plot_it = False):

    if ZIPPERDB_DEBUG:
        print("data_SEG: DEBUG MODE: using CDC19 data")
        # Get the directory of the current script
        current_script_dir = os.path.dirname(__file__)

        # Construct the path to the file
        file_path = os.path.join(current_script_dir, "zipperdb_results_CDC19.csv")
        results_filename = file_path
    else:
            
        # Fetch the results URL
        results_url = fetch_zipperdb_results(protein_fasta)

        # Download the results file

        results_filename = download_results_file(results_url)

    # Parse the results file
    results_df = parse_results_file(results_filename)

    # Create a list to hold the scores
    sequence_only = protein_fasta.seq
    protein_length = len(sequence_only)
    scores = np.zeros(protein_length)

    # Fill the scores list with the scores from the CSV file
    for index, row in results_df.iterrows():
        position = int(row['Position']) - 1  # Convert to 0-based index
        scores[position] = row['Score']

    plot_it = False
    if plot_it:
        # Plotting the results as a heatmap
        plt.figure(figsize=(10, 2))
        plt.imshow(scores.reshape(1, -1), cmap='viridis', aspect='auto')
        plt.colorbar(label='Score')
        plt.yticks([])
        plt.xticks(ticks=np.arange(protein_length), labels=list(sequence_only), fontsize=8)
        plt.xlabel('Amino Acid Position')
        plt.title('Protein Sequence Analysis - CDC19')
        plt.show()
    return scores.reshape(1,-1)


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
    fetch_single_zipperdb_result(cdc19Fasta)
    


if __name__ == "__main__":
    gather_cdc19_data()