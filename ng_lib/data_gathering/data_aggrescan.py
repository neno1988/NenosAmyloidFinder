import requests

import sys
import sys
from pathlib import Path
import numpy as np
from bs4 import BeautifulSoup

# add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ng_lib.data_gathering.interface import DataGatheringTool

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



class AggrescanDataGatheringToolParameters(DataGatheringTool):
    threshold: float = 1e-6

class AggrescanDataGatheringTool(DataGatheringTool):
    def get_data_from_sequence(self, sequence, parameters: AggrescanDataGatheringToolParameters) -> np.ndarray:
        response = self.run_aap(">example\r\n" + sequence)
        print(response)
        result = self.parse_aggrescan_html(response)
        return np.where(np.asarray(result)>parameters.threshold,1,0)


    def run_aap(self, sequence: str):
        url = "http://bioinf.uab.es/cgi-bin/aap/aap_ov.pl"

        data = {
            "sequence": sequence
        }

        # Some servers require a browser-like header to avoid blocking scripts
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
        }

        response = requests.post(url, data=data, headers=headers)

        # Raise error if request failed
        response.raise_for_status()

        return response.text


    def parse_aggrescan_html(self, html: str):

        soup = BeautifulSoup(html, "html.parser")

        def td_spans_to_list(td):
            """Get all span texts from a <td> as a clean list."""
            return [span.get_text(strip=True)
                    for span in td.find_all("span")]

        # 1) Find the header row with '#', 'AA', 'a4v', ...
        header_row = None
        for tr in soup.find_all("tr"):
            cells = [c.get_text(strip=True) for c in tr.find_all(["td", "th"])]
            if cells[:5] == ["#", "AA", "a4v", "HSA", "NHSA"]:
                header_row = tr
                break

        if header_row is None:
            raise RuntimeError("Header row not found")

        # 2) The data row is the next <tr>
        data_row = header_row.find_next_sibling("tr")
        tds = data_row.find_all("td")

        # 3) Extract columns by index
        pos_td   = tds[0]   # positions
        aa_td    = tds[1]   # amino acids
        a4v_td   = tds[2]   # a4v
        hsa_td   = tds[3]   # HSA
        nhsa_td  = tds[4]   # NHSA
        a4vAHS_td = tds[5]  # a4vAHS

        positions = [int(x) for x in td_spans_to_list(pos_td)]
        aa        = td_spans_to_list(aa_td)
        a4v       = [float(x) for x in td_spans_to_list(a4v_td)]
        hsa       = [float(x) for x in td_spans_to_list(hsa_td)]
        nhsa      = [float(x) for x in td_spans_to_list(nhsa_td)]
        a4vAHS    = [float(x) for x in td_spans_to_list(a4vAHS_td)]

        print("positions:", positions)
        print("aa:", aa)
        print("a4v:", a4v)
        return nhsa

if __name__ == "__main__":
    seq = """>example\r\nAMADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERLQSE"""
    tool = AggrescanDataGatheringTool()
    result_html = tool.run_aap(seq)
    result_vector = tool.parse_aggrescan_html(result_html)
    print(result_html)
    print(result_vector)
    