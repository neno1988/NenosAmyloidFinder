import ng_lib as ng
import filecache
import time
import os

import os
import json
import matplotlib.pyplot as plt
import numpy as np
import json
import re
import numpy as np
from ng_lib.visualization.plotly_heatmaps import HeatmapElement, create_interactive_heatmaps

cdc19Fasta = "MSRLERLTSLNVVAGSDLRRTSIIGTIGPKTNNPETLVALRKAGLNIVRMNFSHGSYEYH\
    KSVIDNARKSEELYPGRPLAIALDTKGPEIRTGTTTNDVDYPIPPNHEMIFTTDDKYAKA\
    CDDKIMYVDYKNITKVISAGRIIYVDDGVLSFQVLEVVDDKTLKVKALNAGKICSHKGVN\
    LPGTDVDLPALSEKDKEDLRFGVKNGVHMVFASFIRTANDVLTIREVLGEQGKDVKIIVK\
    IENQQGVNNFDEILKVTDGVMVARGDLGIEIPAPEVLAVQKKLIAKSNLAGKPVICATQM\
    LESMTYNPRPTRAEVSDVGNAILDGADCVMLSGETAKGNYPINAVTTMAETAVIAEQAIA\
    YLPNYDDMRNCTPKPTSTTETVAASAVAAVFEQKAKAIIVLSTSGTTPRLVSKYRPNCPI\
    ILVTRCPRAARFSHLYRGVFPFVFEKEPVSDWTDDVEARINFGIEKAKEFGILKKGDTYV\
    SIQGFKAGAGHSNTLQVSTV"

test_fasta = ">sp|P32119|PRDX2_HUMAN Peroxiredoxin-2 OS=Homo sapiens OX=9606 GN=PRDX2 PE=1 SV=5MASGNARIGKPAPDFKATAVVDGAFKEVKLSDYKGKYVVLFFYPLDFTFVCPTEIIAFSNRAEDFRKLGCEVLGVSVDSQFTHLAWINTPRKEGGLGPLNIPLLADVTRRLSEDYGVLKTDEGIAYRGLFIIDGKGVLRQITVNDLPVGRSVDEALRLVQAFQYTDEHGEVCPAGWKPGSDTIKPNVDDSKEYFSKHN"

# Configuration file path
CONFIG_FILE = 'neno_config.json'

def main_debug():
    
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
    seg_data = np.random.random(len(fasta_sequence.seq))
    zipperDB_data = np.random.random(len(fasta_sequence.seq))
    amylpred_data = np.random.random(len(fasta_sequence.seq))*6


    config = ng.utils.load_config(path = CONFIG_FILE)
    zipperDB_settings = {}
    zipperDB_settings["number_of_rows"] = 2
    zipperDB_settings["number_of_rows"] = 2


    heatmaps: list[HeatmapElement] = []
    heatmaps.append(HeatmapElement(seg_data, number_of_rows=2, min_value=0, max_value=1, colors = ["White", "orange"], legend = "LCR"))
    heatmaps.append(HeatmapElement(amylpred_data, number_of_rows=7, min_value=0, max_value=6, colors = ["White", "green"], legend = "Amylpred"))
    heatmaps.append(HeatmapElement(zipperDB_data, number_of_rows=2, min_value=0, max_value=1, colors = ["White", "darkgray"], legend = "ZipperDB"))
 

    fig = ng.visualization.create_interactive_heatmaps(
        heatmaps=heatmaps,
        name="name",
        xticks=50
    )
    
    ng.visualization.make_annotations(fig, fasta_sequence.name, "left", 25, -23)
    fig.show()




if __name__ == "__main__":
    
    
    config = ng.utils.load_config(path = CONFIG_FILE)
    # Run the application (GUI)
    root = ng.create_gui(config)

    root.mainloop()

    #main_debug() # debugging