import urllib.request
from Bio import SeqIO
from io import StringIO

def get_protein_fasta_by_id(id):
    string_id = "P" + str(id).zfill(5)
    return get_protein_fasta_by_string_id(string_id)

def get_protein_fasta_by_string_id(string_id):
    url = f"https://www.uniprot.org/uniprotkb/{string_id}.fasta"
    return get_protein_fasta(url)

def get_protein_fasta(url):
    # Fetch the content from the URL
    with urllib.request.urlopen(url) as f:
        fasta_data = f.read().decode('utf-8')
    
    # Use StringIO to handle the string as a file object
    fasta_io = StringIO(fasta_data)
    
    # Parse the FASTA format data
    fast_seq = list(SeqIO.parse(fasta_io, "fasta"))
    
    return fast_seq[0]
    

if __name__ == "__main__":
    cdc19_fasta_url = 'https://www.uniprot.org/uniprotkb/P00549.fasta'
    protein_fasta = get_protein_fasta(cdc19_fasta_url)
    print(protein_fasta.seq)
    protein_fasta = get_protein_fasta_by_string_id("P00549")
    print(protein_fasta.seq)
    protein_fasta = get_protein_fasta_by_id(549)
    print(protein_fasta.seq)