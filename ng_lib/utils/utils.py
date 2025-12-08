
import os
import json

CONFIG_FILE = 'neno_config.json'
# Function to save configuration
def save_config(values):
    with open(CONFIG_FILE, 'w') as f:
        json.dump(values, f)

# Function to load configuration
def load_config(path = CONFIG_FILE):
    if os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return {}

# Function to parse FASTA file
def parse_fasta_from_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        description = lines[0].strip()
        sequence = ''.join(line.strip() for line in lines[1:])
        name = description.split(' ')[0][1:]  # Extract ID from description
        return name, description, sequence


class FastaSeq:
    seq:str
    description:str
    name:str

    def __init__(self, seq="", description="", name=""):
        self.seq = seq
        self.description = description
        self.description = name
        