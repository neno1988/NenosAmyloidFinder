import ng_lib as ng
import filecache
import time
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import json
import matplotlib.pyplot as plt
import numpy as np



cdc19Fasta = "MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYH\
    KSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKA\
    CDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVN\
    LPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVK\
    IENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQM\
    LESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIA\
    YLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPI\
    ILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYV\
    SIQGFKAGAGHSNTLQVSTV"


def main():
    somethingFasta= ">sp|P32119|PRDX2_HUMAN Peroxiredoxin-2 OS=Homo sapiens OX=9606 GN=PRDX2 PE=1 SV=5MASGNARIGKPAPDFKATAVVDGAFKEVKLSDYKGKYVVLFFYPLDFTFVCPTEIIAFSNRAEDFRKLGCEVLGVSVDSQFTHLAWINTPRKEGGLGPLNIPLLADVTRRLSEDYGVLKTDEGIAYRGLFIIDGKGVLRQITVNDLPVGRSVDEALRLVQAFQYTDEHGEVCPAGWKPGSDTIKPNVDDSKEYFSKHN"

    Grasp55 = ng.utils.FastaSeq()
    Grasp55.description =  ">sp|Q9H8Y8|GORS2_HUMAN Golgi reassembly-stacking protein 2 OS=Homo sapiens OX=9606 GN=GORASP2 PE=1 SV=3"
    Grasp55.seq = "MGSSQSVEIPGGGTEGYHVLRVQENSPGHRAGLEPFFDFIVSINGSRLNKDNDTLKDLLKANVEKPVKMLIYSSKTLELRETSVTPSNLWGGQGLLGVSIRFCSFDGANENVWHVLEVESNSPAALAGLRPHSDYIIGADTVMNESEDLFSLIETHEAKPLKLYVYNTDTDNCREVIITPNSAWGGEGSLGCGIGYGYLHRIPTRPFEEGKKISLPGQMAGTPITPLKDGFTEVQLSSVNPPSLSPPGTTGIEQSLTGLSISSTPPAVSSVLSTGVPTVPLLPPQVNQSLTSVPPMNPATTLPGLMPLPAGLPNLPNLNLNLPAPHIMPGVGLPELVNPGLPPLPSMPPRNLPGIAPLPLPSEFLPSFPLVPESSSAASSGELLSSLPPTSNAPSDPATTTAKADAASSLTVDVTPPTAKAPTTVEDRVGDSTPVSEKPVSAAVDANASESP"
    Grasp55.name = "Grasp55"

    Grasp65 = ng.utils.FastaSeq()
    Grasp65.description =  ">sp|Q9BQQ3|GORS1_HUMAN Golgi reassembly-stacking protein 1 OS=Homo sapiens OX=9606 GN=GORASP1 PE=1 SV=3"
    Grasp65.seq = "MGLGVSAEQPAGGAEGFHLHGVQENSPAQQAGLEPYFDFIITIGHSRLNKENDTLKALLKANVEKPVKLEVFNMKTMRVREVEVVPSNMWGGQGLLGASVRFCSFRRASEQVWHVLDVEPSSPAALAGLRPYTDYVVGSDQILQESEDFFTLIESHEGKPLKLMVYNSKSDSCREVTVTPNAAWGGEGSLGCGIGYGYLHRIPTQPPSYHKKPPGTPPPSALPLGAPPPDALPPGPTPEDSPSLETGSRQSDYMEALLQAPGSSMEDPLPGPGSPSHSAPDPDGLPHFMETPLQPPPPVQRVMDPGFLDVSGISLLDNSNASVWPSLPSSTELTTTAVSTSGPEDICSSSSSHERGGEATWSGSEFEVSFLDSPGAQAQADHLPQLTLPDSLTSAASPEDGLSAELLEAQAEEEPASTEGLDTGTEAEGLDSQAQISTTE"
    Grasp65.name = "Grasp65"

    fasta_sequence = ng.get_protein_fasta_by_string_id("P32119")
    fasta_sequence = Grasp65
    seg_data = ng.get_SEG_data(fasta_sequence, seg=25)
    zipperDE_data = ng.get_zipperDB_data(fasta_sequence)
    amylpred_data = ng.get_Amylpred_data(fasta_sequence.seq)

    #fig1, _ = ng.nice_heatmap_plot(seg_data)
    #fig1.show()
        # Create interactive plot
    fig = ng.visualization.create_interactive_heatmaps(
        lcr_data=seg_data.reshape(1, -1),
        zipperdb_data=zipperDE_data,
        amylpred_data=amylpred_data if amylpred_data is not None else np.zeros_like(zipperDE_data),
        threshold=-23,
        name="name",
        xticks=50
    )
    fig.show()


@filecache.filecache(24*60*60)
def slow_function(x, y):
    time.sleep(30)
    return x + y
def analyse_protein(output_folder, seq, description, name, threshold, SEG, xticks):
    fasta_protein = ng.utils.FastaSeq()
    fasta_protein.description = description
    fasta_protein.seq = seq
    fasta_protein.name = name

    # Get SEG and ZipperDB data
    seg_data = ng.get_SEG_data(fasta_protein, seg=SEG)
    zipperDE_data = ng.get_zipperDB_data(fasta_protein)

    # Get Amylpred data if credentials are provided
    try:
        amylpred_data = ng.get_Amylpred_data(fasta_protein.seq)
    except Exception as e:
        messagebox.showwarning("Warning", f"Failed to get Amylpred data: {str(e)}")
        amylpred_data = None

    if SEG==0:
        # TODO: only return fig  in nice_heatmap_plot?
        white_cmap = ["white", "white"]
        fig = ng.heatmaps_binary_non_binary(seg_data, zipperDE_data[0], threshold=threshold, name=name, force_cmap=white_cmap, seg_bar_height=0, xticks=xticks)
    else:
        fig = ng.heatmaps_binary_non_binary(seg_data, zipperDE_data[0], threshold=threshold, name=name, xticks=xticks)
    
    
    # Create interactive plot
    fig = ng.visualization.plotly_heatmaps.create_interactive_heatmaps(
        lcr_data=seg_data.reshape(1, -1),
        zipperdb_data=zipperDE_data,
        amylpred_data=amylpred_data if amylpred_data is not None else np.zeros_like(zipperDE_data),
        threshold=threshold,
        name=name,
        xticks=xticks
    )

    # Save the figure as HTML
    html_file = os.path.join(output_folder, f"{name}_analysis.html")
    fig.write_html(html_file)

    # Save the figure as a JPEG file
    fig.savefig(fasta_protein.name+'_treshold_'+str(threshold)+'.jpg', format='jpeg', dpi=300)
    file_name = fasta_protein.name+'_treshold_'+str(threshold)+'.jpg'
    fig.savefig(os.path.join(output_folder, file_name), format='jpeg', dpi=300)
    fig.show()

    # If Amylpred data is available, show it
    if amylpred_data is not None:
        plt.figure(figsize=(10, 2))
        plt.imshow(amylpred_data, cmap='viridis', aspect='auto')
        plt.colorbar(label='Amylpred Score')
        plt.yticks([])
        plt.xticks(ticks=np.arange(len(seq)), fontsize=8)
        plt.xlabel('Amino Acid Position')
        plt.title(f'Amylpred Analysis - {name}')
        plt.show()

### ------------------------------- DATA VALIDATION STUFF ---------------------
import re

def validate_threshold(threshold):
    try:
        threshold = float(threshold)
        if threshold > 0:
            raise ValueError("Threshold must be zero or negative.")
        return threshold
    except ValueError as e:
        raise ValueError(f"Invalid threshold value: {e}")

def validate_sequence(sequence):
    if not re.match("^[ACDEFGHIKLMNPQRSTVWY]*$", sequence):
        raise ValueError("Sequence contains invalid characters. Only amino acid letters are acceptable.")
    return sequence

def validate_xticks(xticks):
    try:
        xticks = int(xticks)
        if xticks < 0:
            raise ValueError("X ticks must be greater or equal to 0.")
        return xticks
    except ValueError as e:
        raise ValueError(f"Invalid X ticks value: {e}")

### ------------------------------- GUI STUFF ---------------------------------------------
# Configuration file path
CONFIG_FILE = 'neno_config.json'

# Function to save configuration
def save_config(values):
    with open(CONFIG_FILE, 'w') as f:
        json.dump(values, f)

# Function to load configuration
def load_config():
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return json.load(f)
    return {}

# Function to parse FASTA file
def parse_fasta(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            description = lines[0].strip()
            sequence = ''.join(line.strip() for line in lines[1:])
            name = description.split(' ')[0][1:]  # Extract ID from description
            return name, description, sequence
    except Exception as e:
        messagebox.showerror("Error", f"Failed to parse FASTA file: {e}")
        return None, None, None

# Function to be called on "Run analysis"
def call_me(output_folder, seq, description, name, threshold, SEG):
    print(f"Output folder: {output_folder}")
    print(f"Sequence: {seq}")
    print(f"Description: {description}")
    print(f"Name: {name}")
    print(f"Threshold: {threshold}")
    print(f"SEG: {SEG}")
    print("Running analysis...")

# Load configuration
config = load_config()

# Create main application window
root = tk.Tk()
root.title("Neno's Amyloid Finder")

# Protein Data Frame
protein_frame = ttk.LabelFrame(root, text="Protein Data")
protein_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

# Name field
name_label = ttk.Label(protein_frame, text="Name:")
name_label.grid(row=0, column=0, sticky="w")
name_entry = ttk.Entry(protein_frame, width=50)
name_entry.grid(row=0, column=1, sticky="ew")
name_entry.insert(0, config.get("name", ""))

# Description field
description_label = ttk.Label(protein_frame, text="Description:")
description_label.grid(row=1, column=0, sticky="w")
description_entry = ttk.Entry(protein_frame, width=50)
description_entry.grid(row=1, column=1, sticky="ew")
description_entry.insert(0, config.get("description", ""))

# Sequence field
sequence_label = ttk.Label(protein_frame, text="Sequence:")
sequence_label.grid(row=2, column=0, sticky="nw")
sequence_text = tk.Text(protein_frame, width=50, height=10, wrap=tk.CHAR)
sequence_text.grid(row=2, column=1, sticky="ew")
sequence_text.insert("1.0", config.get("sequence", ""))

# Scrollbar for Sequence field
sequence_scrollbar = ttk.Scrollbar(protein_frame, orient="vertical", command=sequence_text.yview)
sequence_scrollbar.grid(row=2, column=2, sticky="ns")
sequence_text.configure(yscrollcommand=sequence_scrollbar.set)

# Import from FASTA button
import_button = ttk.Button(protein_frame, text="Import from FASTA", command=lambda: import_fasta())
import_button.grid(row=3, column=1, sticky="e")

# Parameters Frame
parameters_frame = ttk.LabelFrame(root, text="Parameters")
parameters_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

# SEG Dropdown
seg_label = ttk.Label(parameters_frame, text="SEG:")
seg_label.grid(row=0, column=0, sticky="w")
seg_dropdown = ttk.Combobox(parameters_frame, values=["0", "12", "25", "45"], state="readonly")
seg_dropdown.grid(row=0, column=1, sticky="ew")
seg_dropdown.set(config.get("SEG", "0"))

# ZipperDB Threshold field
threshold_label = ttk.Label(parameters_frame, text="ZipperDB Threshold:")
threshold_label.grid(row=1, column=0, sticky="w")
threshold_entry = ttk.Entry(parameters_frame, width=20)
threshold_entry.grid(row=1, column=1, sticky="w")
threshold_entry.insert(0, config.get("threshold", ""))

# Amylpred Credentials Frame
amylpred_frame = ttk.LabelFrame(parameters_frame, text="Amylpred Credentials")
amylpred_frame.grid(row=2, column=0, columnspan=3, padx=5, pady=5, sticky="ew")

# Amylpred Username
amylpred_username_label = ttk.Label(amylpred_frame, text="Username:")
amylpred_username_label.grid(row=0, column=0, sticky="w")
amylpred_username_entry = ttk.Entry(amylpred_frame, width=20)
amylpred_username_entry.grid(row=0, column=1, sticky="w")
amylpred_username_entry.insert(0, config.get("amylpred_username", ""))

# Amylpred Password
amylpred_password_label = ttk.Label(amylpred_frame, text="Password:")
amylpred_password_label.grid(row=1, column=0, sticky="w")
amylpred_password_entry = ttk.Entry(amylpred_frame, width=20, show="*")
amylpred_password_entry.grid(row=1, column=1, sticky="w")
amylpred_password_entry.insert(0, config.get("amylpred_password", ""))

# Output Folder field
output_label = ttk.Label(parameters_frame, text="Output Folder:")
output_label.grid(row=3, column=0, sticky="w")
output_entry = ttk.Entry(parameters_frame, width=40)
output_entry.grid(row=3, column=1, sticky="w")
output_entry.insert(0, config.get("output_folder", ""))
output_button = ttk.Button(parameters_frame, text="Browse", command=lambda: browse_output())
output_button.grid(row=3, column=2, sticky="e")

# Add X ticks parameter
xticks_label = ttk.Label(parameters_frame, text="X ticks (0=auto):").grid(row=4, column=0, padx=5, pady=5, sticky="w")
xticks_var = tk.IntVar(value=0)
xticks_entry = ttk.Entry(parameters_frame, textvariable=xticks_var)
xticks_entry.grid(row=4, column=1, padx=5, pady=5, sticky="ew")

# Run Analysis and Exit Buttons
run_button = ttk.Button(root, text="Run Analysis", command=lambda: run_analysis())
run_button.grid(row=2, column=0, padx=10, pady=5, sticky="ew")
exit_button = ttk.Button(root, text="Exit", command=lambda: exit_app())
exit_button.grid(row=3, column=0, padx=10, pady=5, sticky="ew")

# Helper Functions
def browse_output():
    folder = filedialog.askdirectory()
    if folder:
        output_entry.delete(0, tk.END)
        output_entry.insert(0, folder)

def import_fasta():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta;*.fa")])
    if file_path:
        name, description, sequence = parse_fasta(file_path)
        if name and description and sequence:
            name_entry.delete(0, tk.END)
            name_entry.insert(0, name)
            description_entry.delete(0, tk.END)
            description_entry.insert(0, description)
            sequence_text.delete("1.0", tk.END)
            sequence_text.insert("1.0", sequence)

def run_analysis():
    try:
        threshold = validate_threshold(threshold_entry.get())
        sequence = validate_sequence(sequence_text.get("1.0", tk.END).strip())
        xticks = validate_xticks(xticks_entry.get())
        
        # Set Amylpred credentials in the environment
        os.environ['AMYLPRED_USERNAME'] = amylpred_username_entry.get()
        os.environ['AMYLPRED_PASSWORD'] = amylpred_password_entry.get()
   
        analyse_protein(
            output_folder=output_entry.get(),
            seq=sequence,
            description=description_entry.get(),
            name=name_entry.get(),
            threshold=threshold,
            SEG=int(seg_dropdown.get()), 
            xticks=xticks
        )
        messagebox.showinfo("Success", "Analysis complete!")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to run analysis: {e}")

def exit_app():
    save_config({
        "name": name_entry.get(),
        "description": description_entry.get(),
        "sequence": sequence_text.get("1.0", tk.END).strip(),
        "SEG": seg_dropdown.get(),
        "threshold": threshold_entry.get(),
        "output_folder": output_entry.get(),
        "xticks": xticks_entry.get(),
        "amylpred_username": amylpred_username_entry.get(),
        "amylpred_password": amylpred_password_entry.get()
    })
    root.destroy()



if __name__ == "__main__":
    # Run the application (GUI)
    #root.mainloop()

    main() # debugging