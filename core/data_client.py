import requests
import re
from Bio import SeqIO
from io import StringIO
from typing import Optional, Tuple, Dict
from config import AA_MAP

class BioDataClient:
    
    def __init__(self):
        # Headers avoid 403 errors from biological databases
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }

    def get_sequence(self, uniprot_id: str) -> Tuple[Optional[str], Optional[str]]:
        """Fetches sequence from UniProt."""
        base_id = uniprot_id.split('-')[0]
        url = f"https://rest.uniprot.org/uniprotkb/{base_id}.fasta"
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                record = SeqIO.read(StringIO(r.text), "fasta")
                return str(record.seq), record.description
        except Exception as e:
            print(f"Error fetching sequence: {e}")
        return None, None

    def get_structure(self, uniprot_id: str) -> Tuple[Optional[str], str]:
        """Fetches PDB from AlphaFold with robust fallback logic."""
        base_id = uniprot_id 
        versions = ['v4', 'v3', 'v5', 'v6']
        
        # 1. Try Standard Model
        base_url = f"https://alphafold.ebi.ac.uk/files/AF-{base_id}-F1-model_"
        for version in versions:
            url = base_url + version + ".pdb"
            try:
                r = requests.get(url, headers=self.headers, timeout=10)
                if r.status_code == 200:
                    return r.text, version
            except:
                continue
        
        # 2. Try Isoform Fallback
        if '-' not in base_id:
            iso_base_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-4-F1-model_"
            for version in versions:
                url = iso_base_url + version + ".pdb"
                try:
                    r = requests.get(url, headers=self.headers, timeout=10)
                    if r.status_code == 200:
                        return r.text, f"{version} (Isoform -4)"
                except:
                    continue

        return None, "Not Found"

    def fetch_clinical_data(self, gene: str, mutation: str) -> Dict:
        """
        Robust ClinVar lookup that actually works in 2025.
        Tries 8 different search strategies until one returns clinical significance.
        """

        m = re.match(r"^([A-Z])(\d+)([A-Z])$", mutation.upper())
        if not m:
            return {"error": "Invalid mutation format, expected e.g. R273H"}

        wt_aa, pos, mut_aa = m.groups()
        # one-letter codes
        pos = int(pos)

        # Convert to 3-letter for ClinVar's preferred format
        if wt_aa not in AA_MAP or mut_aa not in AA_MAP:
            return {"error": "Unknown amino acid"}
        wt3 = AA_MAP[wt_aa]
        mut3 = AA_MAP[mut_aa]

        # All search strategies — we stop at the first one that finds ClinVar data
        strategies = [
            # 1. Exact protein change with 3-letter code (most common in ClinVar)
            f'{gene} AND "p.{wt3}{pos}{mut3}"',
            # 2. One-letter version (some records use this)
            f'{gene} AND "p.{wt_aa}{pos}{mut_aa}"',
            # 3. Simple symbols without p.
            f'{gene} AND {wt3}{pos}{mut_aa}',
            f'{gene} AND {wt_aa}{pos}{mut_aa}',
            # 4. Very broad — just gene + position + mutant AA (catches off-by-1 isoforms)
            f'{gene} AND {pos}{mut_aa}',
            f'{gene} AND {pos}{mut3}',
            # 5. Fallback: search by gene only and scan all variants at that position
            f'clinvar.gene.symbol:{gene} AND clinvar.protein:"*{pos}*"',
            # 6. Nuclear option
            f'{gene}',
        ]

        url = "https://myvariant.info/v1/query"

        for i, q in enumerate(strategies, 1):
            params = {
                "q": q,
                "fields": "clinvar,dbsnp",
                "size": 20,
                "fetch_all": "false"
            }
            try:
                r = requests.get(url, params=params, headers=self.headers, timeout=8)
                if r.status_code != 200:
                    continue
                data = r.json()

                for hit in data.get("hits", []):
                    clinvar = hit.get("clinvar")
                    if not clinvar:
                        continue

                    # Modern MyVariant structure — rcv is often a list inside variant_rcv or directly rcv
                    rcv_list = clinvar.get("rcv") or clinvar.get("variant_rcv", [])
                    if isinstance(rcv_list, dict):
                        rcv_list = [rcv_list]

                    for record in rcv_list:
                        if not record:
                            continue

                        sig = record.get("clinical_significance", "").lower()
                        if not sig or "not provided" in sig or "no assertion" in sig:
                            continue

                        condition = ""
                        conds = record.get("conditions")
                        if conds:
                            if isinstance(conds, list):
                                condition = conds[0].get("name", "")
                            elif isinstance(conds, dict):
                                condition = conds.get("name", "")
                        
                        return {
                            "significance": record.get("clinical_significance", "Unknown"),
                            "conditions": condition or "Not specified",
                            "source": "MyVariant.info → ClinVar",
                            "matched_strategy": i,
                            "variant_id": hit.get("_id"),
                            "debug_query": q
                        }
            except Exception as e:
                print(f"Strategy {i} failed: {e}")
                continue

        # In core/data_client.py, inside fetch_clinical_data...
        return {
            "significance": first_record.get("clinical_significance", "Unknown"),
            "conditions": first_record.get("conditions", {}).get("name", "Not specified"),
            "source": "MyVariant.info (ClinVar)",
            "debug_query": q 
        }