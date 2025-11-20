# config.py
from typing import Dict

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"

KNOWN_MUTATIONS: Dict[str, Dict[str, str]] = {
    'TP53_R175H': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'Hotspot mutation inactivating p53 tumor suppressor function.'
    },
    'TP53_R248Q': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'Disrupts DNA-binding domain, loss of transcriptional activity.'
    },
    'TP53_R273H': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'Gain-of-function properties, enhancing cancer cell survival.'
    }
}

# Standard Amino Acid Mapping
AA_MAP = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}