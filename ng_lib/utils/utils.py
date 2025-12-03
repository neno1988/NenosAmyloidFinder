
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


class FastaSeq:
    seq:str
    description:str
    name:str

    def __init__(self, seq="", description="", name=""):
        self.seq = seq
        self.description = description
        self.description = name
        