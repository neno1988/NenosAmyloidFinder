import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import json
import os
from ng_lib.utils import FastaSeq, save_config, load_config, parse_fasta_from_file
from ng_lib.analysis import analyse_protein

def create_gui(config):
    # Create main application window
    root = tk.Tk()
    root.title("Neno's Amyloid Finder")

    # Protein Data Frame
    root.protein_frame = ttk.LabelFrame(root, text="Protein Data")
    root.protein_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

    # Name field
    root.name_label = ttk.Label(root.protein_frame, text="Name:")
    root.name_label.grid(row=0, column=0, sticky="w")
    root.name_entry = ttk.Entry(root.protein_frame, width=50)
    root.name_entry.grid(row=0, column=1, sticky="ew")
    root.name_entry.insert(0, config.get("name", ""))
    # Description field
    root.description_label = ttk.Label(root.protein_frame, text="Description:")
    root.description_label.grid(row=1, column=0, sticky="w")
    root.description_entry = ttk.Entry(root.protein_frame, width=50)
    root.description_entry.grid(row=1, column=1, sticky="ew")
    root.description_entry.insert(0, config.get("description", ""))

    # Sequence field
    root.sequence_label = ttk.Label(root.protein_frame, text="Sequence:")
    root.sequence_label.grid(row=2, column=0, sticky="nw")
    root.sequence_text = tk.Text(root.protein_frame, width=50, height=10, wrap=tk.CHAR)
    root.sequence_text.grid(row=2, column=1, sticky="ew")
    root.sequence_text.insert("1.0", config.get("sequence", ""))

    # Scrollbar for Sequence field
    root.sequence_scrollbar = ttk.Scrollbar(root.protein_frame, orient="vertical", command=root.sequence_text.yview)
    root.sequence_scrollbar.grid(row=2, column=2, sticky="ns")
    root.sequence_text.configure(yscrollcommand=root.sequence_scrollbar.set)

    # Import from FASTA button
    root.import_button = ttk.Button(root.protein_frame, text="Import from FASTA", command=lambda: import_fasta(root))
    root.import_button.grid(row=3, column=1, sticky="e")
    # Parameters Frame
    root.parameters_frame = ttk.LabelFrame(root, text="Parameters")
    root.parameters_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

    # SEG Dropdown
    root.seg_label = ttk.Label(root.parameters_frame, text="SEG:")
    root.seg_label.grid(row=0, column=0, sticky="w")
    root.seg_dropdown = ttk.Combobox(root.parameters_frame, values=["0", "12", "25", "45"], state="readonly")
    root.seg_dropdown.grid(row=0, column=1, sticky="ew")
    root.seg_dropdown.set(config.get("SEG", "0"))

    # ZipperDB Threshold field
    root.threshold_label = ttk.Label(root.parameters_frame, text="ZipperDB Threshold:")
    root.threshold_label.grid(row=1, column=0, sticky="w")
    root.threshold_entry = ttk.Entry(root.parameters_frame, width=20)
    root.threshold_entry.grid(row=1, column=1, sticky="w")
    root.threshold_entry.insert(0, config.get("threshold", ""))

    # Amylpred Credentials Frame
    root.amylpred_frame = ttk.LabelFrame(root.parameters_frame, text="Amylpred Credentials")
    root.amylpred_frame.grid(row=2, column=0, columnspan=3, padx=5, pady=5, sticky="ew")
    # Amylpred Username
    root.amylpred_username_label = ttk.Label(root.amylpred_frame, text="Username:")
    root.amylpred_username_label.grid(row=0, column=0, sticky="w")
    root.amylpred_username_entry = ttk.Entry(root.amylpred_frame, width=20)
    root.amylpred_username_entry.grid(row=0, column=1, sticky="w")
    root.amylpred_username_entry.insert(0, config.get("amylpred_username", ""))

    # Amylpred Password
    root.amylpred_password_label = ttk.Label(root.amylpred_frame, text="Password:")
    root.amylpred_password_label.grid(row=1, column=0, sticky="w")
    root.amylpred_password_entry = ttk.Entry(root.amylpred_frame, width=20, show="*")
    root.amylpred_password_entry.grid(row=1, column=1, sticky="w")
    root.amylpred_password_entry.insert(0, config.get("amylpred_password", ""))

    # Output Folder field
    root.output_label = ttk.Label(root.parameters_frame, text="Output Folder:")
    root.output_label.grid(row=3, column=0, sticky="w")
    root.output_entry = ttk.Entry(root.parameters_frame, width=40)
    root.output_entry.grid(row=3, column=1, sticky="w")
    root.output_entry.insert(0, config.get("output_folder", ""))
    root.output_button = ttk.Button(root.parameters_frame, text="Browse", command=lambda: browse_output(root))
    root.output_button.grid(row=3, column=2, sticky="e")
    # Add X ticks parameter
    root.xticks_label = ttk.Label(root.parameters_frame, text="X ticks (0=auto):")
    root.xticks_label.grid(row=4, column=0, padx=5, pady=5, sticky="w")
    root.xticks_var = tk.IntVar(value=0)
    root.xticks_entry = ttk.Entry(root.parameters_frame, textvariable=root.xticks_var)
    root.xticks_entry.grid(row=4, column=1, padx=5, pady=5, sticky="ew")

    # Run Analysis and Exit Buttons
    root.run_button = ttk.Button(root, text="Run Analysis", command=lambda: run_analysis(root))
    root.run_button.grid(row=2, column=0, padx=10, pady=5, sticky="ew")
    root.exit_button = ttk.Button(root, text="Exit", command=lambda: exit_app(root))
    root.exit_button.grid(row=3, column=0, padx=10, pady=5, sticky="ew")

    return root

# Helper Functions
def browse_output(root):
    folder = filedialog.askdirectory()
    if folder:
        root.output_entry.delete(0, tk.END)
        root.output_entry.insert(0, folder)

def import_fasta(root):
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta;*.fa")])
    if file_path:
        name, description, sequence = gui_parse_fasta(file_path)
        if name and description and sequence:
            root.name_entry.delete(0, tk.END)
            root.name_entry.insert(0, name)
            root.description_entry.delete(0, tk.END)
            root.description_entry.insert(0, description)
            root.sequence_text.delete("1.0", tk.END)
            root.sequence_text.insert("1.0", sequence)
            
def parse_fasta(fasta_str):
    lines = fasta_str.split("\n")
    description = lines[0].strip()
    sequence = ''.join(line.strip() for line in lines[1:])
    name = description.split(' ')[0][0:]  # Extract ID from description
    ret = FastaSeq()
    ret.name = name
    ret.seq = sequence
    ret.description = description
    return ret

def run_analysis(root):
    #try:

    threshold = validate_threshold(root.threshold_entry.get())
    
    sequence_txt_stripped = root.sequence_text.get("1.0", tk.END).strip()
    if sequence_txt_stripped.count('>')==1:
        # if is fasta format (TODO: use some better function in biotools to detect fasta?)
        fasta_seq = parse_fasta(sequence_txt_stripped)
        sequence = fasta_seq.seq
        sequences = [sequence]
        names = [root.name_entry.get()]
    elif sequence_txt_stripped.count('>')==0:
        # only sequence
        sequence = sequence_txt_stripped.replace("\n", "").replace(" ", "")
        sequences = [sequence]
        names = [root.name_entry.get()]
    elif sequence_txt_stripped.count('>') > 1:
        fastas = sequence_txt_stripped.split('>')
        sequences = []
        names = []
        for f in fastas:
            fasta_object = parse_fasta(f)
            if len(fasta_object.seq) >2:
                sequences.append(fasta_object.seq)
                names.append(fasta_object.name)

    for sequence, name in zip(sequences, names):
        sequence = validate_sequence(sequence)

        xticks = validate_xticks(root.xticks_entry.get())
        
        # Set Amylpred credentials in the environment
        os.environ['AMYLPRED_USERNAME'] = root.amylpred_username_entry.get()
        os.environ['AMYLPRED_PASSWORD'] = root.amylpred_password_entry.get()
        analyse_protein(
            output_folder=root.output_entry.get(),
            seq=sequence,
            description=root.description_entry.get(),
            name=name,
            threshold=threshold,
            SEG=int(root.seg_dropdown.get()), 
            xticks=xticks
        )
    messagebox.showinfo("Success", "Analysis complete!")
    #except Exception as e:
        #messagebox.showerror("Error", f"Failed to run analysis: {e}")

# Function to parse FASTA file
def gui_parse_fasta(file_path):
    try:
        name, description, sequence =parse_fasta_from_file(file_path)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to parse FASTA file: {e}")
        name, description, sequence =  None, None, None
    return name, description, sequence
    
def exit_app(root):
    save_config({
        "name": root.name_entry.get(),
        "description": root.description_entry.get(),
        "sequence": root.sequence_text.get("1.0", tk.END).strip(),
        "SEG": root.seg_dropdown.get(),
        "threshold": root.threshold_entry.get(),
        "output_folder": root.output_entry.get(),
        "xticks": root.xticks_entry.get(),
        "amylpred_username": root.amylpred_username_entry.get(),
        "amylpred_password": root.amylpred_password_entry.get()
    })
    root.destroy()



### ------------------------------- DATA VALIDATION  ---------------------
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