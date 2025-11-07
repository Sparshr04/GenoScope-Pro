import streamlit as st
import requests
from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch
from Bio import SeqIO
from io import StringIO
import re
import py3Dmol
from stmol import showmol
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Optional

st.set_page_config(page_title="GenoScope", page_icon="🧬", layout="wide")

# ----------- MODEL SETUP -----------
@st.cache_resource
def load_model():
    model_name = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    return tokenizer, model

tokenizer, model = load_model()

# ----------- KNOWN MUTATIONS DICT FOR SEVERITY AND SUMMARY -----------
known_mutations: Dict[str, Dict[str, str]] = {
    'R175H': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'The R175H mutation in TP53 is a hotspot mutation that inactivates p53\'s tumor suppressor function and is associated with Li-Fraumeni syndrome and various cancers.'
    },
    'R248Q': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'R248Q disrupts the DNA-binding domain of p53, leading to loss of transcriptional activity and increased cancer risk.'
    },
    'R249S': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'This mutation affects the DNA contact residue in p53, impairing its ability to bind DNA and suppress tumors.'
    },
    'R273H': {
        'severity': '🔴 Likely Deleterious',
        'summary': 'The R273H mutation in p53 confers gain-of-function properties, enhancing cancer cell survival, anoikis resistance, and drug resistance, associated with multiple cancers.'
    },
    # Add more as needed
}

# ----------- FUNCTIONS -----------
def get_sequence(uniprot_id: str) -> tuple[Optional[str], Optional[str]]:
    """Fetch FASTA sequence and description from UniProt"""
    base_id = uniprot_id.split('-')[0]
    url = f"https://rest.uniprot.org/uniprotkb/{base_id}.fasta"
    r = requests.get(url)
    if r.status_code != 200:
        return None, None
    record = SeqIO.read(StringIO(r.text), "fasta")
    return str(record.seq), record.description

def get_uniprot_info(uniprot_id: str) -> Optional[Dict]:
    """Fetch additional UniProt info like function, length"""
    base_id = uniprot_id.split('-')[0]
    url = f"https://rest.uniprot.org/uniprotkb/{base_id}.json"
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        return {
            'length': data['sequence']['length'],
            'function': data.get('comments', [{}])[0].get('texts', [{}])[0].get('value', 'No function description') if data.get('comments') else 'No function description',
            'gene_name': data['genes'][0]['geneName']['value'] if data.get('genes') else 'Unknown'
        }
    return None

def compute_mutation_score(seq: str, mutation: str) -> Optional[float]:
    m = re.match(r"([A-Z])(\d+)([A-Z])", mutation.upper())
    if not m:
        st.error("❌ Mutation must be like A23T")
        return None
    wt, pos_str, mut = m.groups()
    pos = int(pos_str) - 1
    if pos < 0 or pos >= len(seq):
        st.error("❌ Invalid position for this protein.")
        return None
    masked = seq[:pos] + "<mask>" + seq[pos+1:]
    inputs = tokenizer(masked, return_tensors="pt")
    with torch.no_grad():
        logits = model(**inputs).logits
    mask_token_id = tokenizer.mask_token_id
    mask_indices = (inputs.input_ids == mask_token_id).nonzero(as_tuple=True)[1]
    if len(mask_indices) == 0:
        st.error("❌ Mask token not found in input.")
        return None
    mask_pos = mask_indices[0].item()
    probs = torch.nn.functional.softmax(logits[0, mask_pos, :], dim=-1)
    wt_token_id = tokenizer.convert_tokens_to_ids(wt)
    mut_token_id = tokenizer.convert_tokens_to_ids(mut)
    if wt_token_id is None or mut_token_id is None:
        st.error("❌ Invalid amino acid tokens.")
        return None
    wt_score = probs[wt_token_id].item()
    mut_score = probs[mut_token_id].item()
    if wt_score == 0:
        st.warning("⚠️ Wild-type score is zero; delta may be unreliable.")
        return None
    delta = float(torch.log(torch.tensor(mut_score / wt_score)))
    return delta

def classify_impact(delta: float, mutation: str, gene: str = 'TP53') -> str:
    """Classify impact, prioritizing known mutations"""
    key = f"{gene}_{mutation.upper()}"
    if key in known_mutations:
        return known_mutations[key]['severity']
    if delta < -0.5:
        return "🔴 Likely Deleterious"
    elif delta < 0:
        return "🟠 Possibly Harmful"
    else:
        return "🟢 Likely Benign"

def get_mutation_insight(mutation: str, gene: str = 'TP53') -> str:
    """Get AI-powered summary, using known dict for now"""
    key = f"{gene}_{mutation.upper()}"
    if key in known_mutations:
        return known_mutations[key]['summary']
    return f"The {mutation} mutation in {gene} may alter protein function. Further analysis recommended."

# Mock ClinVar fetch - in production, use proper API
def get_clinvar_info(mutation: str, gene: str = 'TP53') -> Optional[Dict]:
    """Mock ClinVar data"""
    mock_data = {
        'R273H': {'clinical_significance': 'Pathogenic', 'conditions': 'Li-Fraumeni syndrome, various cancers', 'rs_id': 'rs28934574'},
        # Add more
    }
    key = mutation.upper()
    if key in mock_data:
        return mock_data[key]
    return None

def get_pdb_from_alphafold(uniprot_id: str) -> tuple[Optional[str], Optional[str]]:
    """Fetch PDB model directly from AlphaFold, trying different versions and isoforms"""
    if '-' in uniprot_id:
        base_id = uniprot_id
    else:
        base_id = uniprot_id
    base_url = f"https://alphafold.ebi.ac.uk/files/AF-{base_id}-F1-model_"
    versions = ['v4', 'v3', 'v5', 'v6']
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    for version in versions:
        url = base_url + version + ".pdb"
        r = requests.get(url, headers=headers)
        if r.status_code == 200:
            st.info(f"✅ Fetched model version {version} for {base_id}")
            return r.text, version
        else:
            st.warning(f"⚠️ Version {version} not found for {base_id} (status: {r.status_code})")
    if '-' not in base_id:
        iso_base_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-4-F1-model_"
        for version in versions:
            url = iso_base_url + version + ".pdb"
            r = requests.get(url, headers=headers)
            if r.status_code == 200:
                st.info(f"✅ Fetched isoform -4 model version {version}")
                return r.text, version
            else:
                st.warning(f"⚠️ Isoform -4 version {version} not found (status: {r.status_code})")
    return None, None

def render_structure(pdb_data: str, mutation: Optional[str] = None):
    """Render PDB in 3Dmol using stmol with mutation highlight"""
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    if mutation:
        m = re.match(r"[A-Z](\d+)[A-Z]", mutation.upper())
        if m:
            pos = int(m.group(1))
            view.addStyle({"resi": str(pos)}, {"stick": {"color": "red", "radius": 0.3}})
            view.zoomTo({"resi": str(pos)})
    showmol(view, height=600, width=800)

def extract_plddt_info(pdb_data: str) -> Optional[Dict]:
    """Parse AlphaFold PDB for pLDDT scores (stored in B-factor column) with residue numbers"""
    if not pdb_data:
        return None
    lines = pdb_data.splitlines()
    res_plddt = {}
    for line in lines:
        if line.startswith("ATOM") and len(line) > 66:
            if line[12:16].strip() == "CA":
                resnum = line[22:26].strip()
                try:
                    plddt = float(line[60:66])
                    if 0 <= plddt <= 100:
                        res_plddt[resnum] = plddt
                except ValueError:
                    pass
    if res_plddt:
        resnums = sorted(int(k) for k in res_plddt.keys())
        plddts = [res_plddt[str(r)] for r in resnums]
        avg_plddt = sum(plddts) / len(plddts)
        df = pd.DataFrame({'Residue': resnums, 'pLDDT': plddts})
        return {
            "df": df,
            "average_plddt": round(avg_plddt, 2),
            "confidence_summary": "🟢 Very High (>90)" if avg_plddt > 90 else "🟡 High (70-90)" if avg_plddt > 70 else "🟠 Low (50-70)" if avg_plddt > 50 else "🔴 Very Low (<50)"
        }
    return None

def analyze_batch_mutations(uniprot_id: str, mutations_file):
    """Analyze batch of mutations"""
    seq, _ = get_sequence(uniprot_id)
    if not seq:
        return
    mutations_text = mutations_file.read().decode("utf-8")
    mutations = [m.strip() for m in mutations_text.splitlines() if m.strip()]
    results = []
    for mut in mutations:
        delta = compute_mutation_score(seq, mut)
        impact = classify_impact(delta if delta is not None else 0, mut)
        insight = get_mutation_insight(mut)
        results.append({'Mutation': mut, 'Impact': impact, 'Delta': delta, 'Insight': insight})
    df = pd.DataFrame(results)
    st.dataframe(df)
    # Simple chart: severity count
    severity_counts = df['Impact'].value_counts()
    st.bar_chart(severity_counts)

# ----------- STREAMLIT UI -----------
st.title("🧬 GenoScope – AI Mutation Impact Predictor")
st.markdown(
    "### Analyze how a genetic mutation affects protein function using **AI-powered prediction** and **3D visualization**."
)
st.divider()

# Input
col1, col2 = st.columns([2, 1])

with col1:
    uniprot_id = st.text_input("Enter UniProt ID", placeholder="e.g. P04637 (TP53)")
    mutation = st.text_input("Enter Mutation", placeholder="e.g. R273H")
    gene = st.text_input("Gene (optional)", value="TP53", placeholder="e.g. TP53")

with col2:
    st.markdown("<div style='display: flex; justify-content: center; align-items: center; height: 100%;'>", unsafe_allow_html=True)
    st.image(
        "https://upload.wikimedia.org/wikipedia/commons/0/0c/DNA_animation.gif",
        width=180,  # This is ~50% of original (~240px) — adjust as needed
        caption="AI meets Genomics"
    )
    st.markdown("</div>", unsafe_allow_html=True)
# Batch upload
uploaded_file = st.file_uploader("📁 Upload batch mutations (one per line)", type=['txt'])
if uploaded_file:
    if st.button("🚀 Analyze Batch"):
        analyze_batch_mutations(uniprot_id, uploaded_file)

if st.button("🚀 Analyze Single Mutation"):
    with st.spinner("Analyzing mutation..."):
        seq, protein_name = get_sequence(uniprot_id)
        uniprot_info = get_uniprot_info(uniprot_id)
        if seq:
            st.session_state['seq'] = seq
            st.session_state['protein_name'] = protein_name
            st.session_state['uniprot_info'] = uniprot_info
            delta = compute_mutation_score(seq, mutation)
            if delta is not None:
                impact = classify_impact(delta, mutation, gene)
                insight = get_mutation_insight(mutation, gene)
                clinvar_info = get_clinvar_info(mutation, gene)
                st.session_state['delta'] = delta
                st.session_state['impact'] = impact
                st.session_state['insight'] = insight
                st.session_state['clinvar'] = clinvar_info
            else:
                st.error("Unable to compute score.")
        else:
            st.error("Invalid UniProt ID.")

# Tabs
tab1, tab2, tab3, tab4 = st.tabs(["Overview", "3D Structure", "Mutation Impact", "Database Info"])

with tab1:
    if 'seq' in st.session_state:
        st.markdown("### 🧬 Protein Sequence")
        st.text_area("Sequence:", value=st.session_state['seq'], height=100, disabled=True)
        if st.session_state.get('uniprot_info'):
            info = st.session_state['uniprot_info']
            col_a, col_b, col_c = st.columns(3)
            with col_a:
                st.metric("Length", info['length'])
            with col_b:
                st.caption("**Function:**")
                st.write(info['function'][:200] + "...")
            with col_c:
                st.caption("**Gene:**")
                st.write(info['gene_name'])
        st.markdown(f"**Protein Name:** {st.session_state['protein_name']}")
        if 'impact' in st.session_state:
            st.success(f"**Predicted Impact:** {st.session_state['impact']}")
            st.info(f"Δ log-probability = {st.session_state['delta']:.3f}")

with tab2:
    if uniprot_id and st.button("🧩 Load 3D Structure", key="load_struct"):
        with st.spinner("Fetching AlphaFold model..."):
            pdb_data, version = get_pdb_from_alphafold(uniprot_id)
            if pdb_data:
                st.session_state['pdb_data'] = pdb_data
                st.session_state['version'] = version
                st.markdown("### 🧠 Protein Structure Visualization (AlphaFold Model)")
                col_left, col_right = st.columns([3, 1])
                with col_left:
                    render_structure(pdb_data, mutation)
                with col_right:
                    st.markdown("#### Model Confidence")
                    st.markdown("**pLDDT (per-residue measure of local confidence)**")
                    legend_data = {
                        "Very high (pLDDT >90)": "🟦",
                        "High (90 > pLDDT >70)": "🔷",
                        "Low (70 > pLDDT >50)": "🟨",
                        "Very low (pLDDT <50)": "🟠"
                    }
                    for label, color in legend_data.items():
                        st.markdown(f"{color} {label}")
                    st.markdown("*Learn more...*")
                    st.markdown("#### Domains (0)")
                    st.markdown("#### Annotations")
                    model_label = f"AF-{uniprot_id}-F1 Model | Instance 1/1 | {version}"
                    st.markdown(f"**{model_label}**")
                
                plddt_info = extract_plddt_info(pdb_data)
                if plddt_info:
                    st.markdown("### 📊 pLDDT Confidence Graph")
                    fig, ax = plt.subplots(figsize=(10, 4))
                    ax.plot(plddt_info["df"]["Residue"], plddt_info["df"]["pLDDT"], color='blue', linewidth=1)
                    ax.set_xlabel("Residue Number")
                    ax.set_ylabel("pLDDT Score")
                    ax.set_title("pLDDT per Residue")
                    ax.fill_between(plddt_info["df"]["Residue"], plddt_info["df"]["pLDDT"], alpha=0.3)
                    ax.axhline(y=90, color='green', linestyle='--', label='Very High')
                    ax.axhline(y=70, color='orange', linestyle='--', label='High')
                    ax.axhline(y=50, color='red', linestyle='--', label='Low')
                    ax.legend()
                    st.pyplot(fig)
                    
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("Average pLDDT", plddt_info["average_plddt"])
                    with col_b:
                        st.info(plddt_info["confidence_summary"])
            else:
                st.error("Structure not found.")
    else:
        st.info("👆 Enter UniProt ID and analyze to load structure.")

with tab3:
    if 'impact' in st.session_state:
        st.markdown("### 🎯 Mutation Impact")
        st.success(st.session_state['impact'])
        st.markdown("### 🧠 AI-Powered Insight")
        st.info(st.session_state['insight'])
    else:
        st.info("👆 Analyze a mutation first.")

with tab4:
    if 'uniprot_info' in st.session_state:
        st.markdown("### 🔗 UniProt Information")
        info = st.session_state['uniprot_info']
        st.json(info)
    if 'clinvar' in st.session_state and st.session_state['clinvar']:
        st.markdown("### 🏥 ClinVar Annotation")
        st.json(st.session_state['clinvar'])
    else:
        st.info("👆 Analyze a mutation for database info.")

st.info("💡 **Batch Mode:** Upload a TXT file with one mutation per line (e.g., R273H\nR175H) for comparative analysis.")