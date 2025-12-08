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
    seg_data = ng.get_SEG_data(fasta_sequence, seg=25)
    zipperDE_data = ng.get_zipperDB_data(fasta_sequence)
    amylpred_data = ng.get_amylpred_data(fasta_sequence.seq)

    #fig1, _ = ng.nice_heatmap_plot(seg_data)
    #fig1.show()
        # Create interactive plot
    debugging_ng = False
    if debugging_ng:
        fig = ng.visualization.create_interactive_heatmaps(
            lcr_data=seg_data.reshape(1, -1),
            seg_treshold = 25,
            zipperdb_data=zipperDE_data,
            amylpred_data=amylpred_data if amylpred_data is not None else np.zeros_like(zipperDE_data),
            zipperdb_threshold=-23,
            name="name",
            xticks=50
        )
        fig.show()
    else: 
        # Create interactive plot
        fig = ng.visualization.plotly_heatmaps.create_interactive_heatmaps(
            lcr_data=seg_data.reshape(1, -1),
            seg_treshold = 25,
            zipperdb_data=zipperDE_data,
            amylpred_data=amylpred_data if amylpred_data is not None else np.zeros_like(zipperDE_data),
            zipperdb_threshold=-23,
            name="name",
            xticks=0
        )


if __name__ == "__main__":
    # Configuration file path
    CONFIG_FILE = 'neno_config.json'
    config = ng.utils.load_config(path = CONFIG_FILE)
    # Run the application (GUI)
    root = ng.create_gui(config)

    root.mainloop()

    #main_debug() # debugging