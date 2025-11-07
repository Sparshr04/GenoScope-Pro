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
import time
import streamlit.components.v1 as components

# ========================= PAGE CONFIG =========================
st.set_page_config(
    page_title="GenoScope • AI Mutation Analyzer",
    page_icon="DNA",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ========================= CLEAN, SAFE CSS (NO RAW DIVS) =========================
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap');
    
    * { font-family: 'Inter', sans-serif; }
    
    .main { 
        background: linear-gradient(135deg, #141e30, #243b55);
        background-attachment: fixed;
    }
    
    h1, h2, h3 { 
        background: linear-gradient(90deg, #a8edea, #fed6e3);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-weight: 700;
    }
    
    .stButton > button {
        background: linear-gradient(45deg, #00c9ff, #92fe9d) !important;
        color: black !important;
        border-radius: 50px !important;
        padding: 12px 32px !important;
        font-weight: 600 !important;
        box-shadow: 0 8px 20px rgba(0,201,255,0.4) !important;
    }
    
    .css-1v0mbdj {  /* Input fields */
        background: rgba(255,255,255,0.1) !important;
        border-radius: 12px !important;
        border: 1px solid rgba(255,255,255,0.2) !important;
    }
    
    .metric-card {
        background: rgba(255,255,255,0.15);
        border-radius: 16px;
        padding: 1.5rem;
        text-align: center;
        backdrop-filter: blur(10px);
        border: 1px solid rgba(255,255,255,0.2);
        box-shadow: 0 8px 32px rgba(0,0,0,0.3);
    }
    
    .insight-box {
        background: rgba(0, 201, 255, 0.1);
        border-left: 5px solid #00c9ff;
        padding: 1rem;
        border-radius: 0 12px 12px 0;
    }
</style>
""", unsafe_allow_html=True)

# ========================= CONFETTI =========================
def celebrate():
    components.html("""
    <script src="https://cdn.jsdelivr.net/npm/canvas-confetti@1.6.0/dist/confetti.browser.min.js"></script>
    <script>
        setTimeout(() => {
            confetti({ particleCount: 150, spread: 70, origin: { y: 0.6 } });
        }, 500);
    </script>
    """, height=0)

# ========================= MODEL =========================
@st.cache_resource
def load_model():
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    model = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
    return tokenizer, model

tokenizer, model = load_model()

# ========================= DATA =========================
INSIGHTS = {
    "R273H": "Gain-of-function hotspot in p53. Disrupts DNA binding and confers dominant-negative effect. Strongly associated with Li-Fraumeni syndrome and multiple cancers.",
    "R175H": "Zinc-binding residue mutation. Destabilizes p53 core domain. Common in ovarian and lung cancers.",
    "R248Q": "Direct DNA contact site. Loss of transcriptional activation. High prevalence in breast and colorectal cancers.",
    "R249S": "Aflatoxin B1-induced mutation. Liver cancer hotspot in high-exposure regions."
}

# ========================= FUNCTIONS =========================
def get_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id.split('-')[0]}.fasta"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            record = SeqIO.read(StringIO(r.text), "fasta")
            return str(record.seq), record.description.split("OS=")[0].strip()
    except:
        pass
    return None, None

def compute_delta(seq, mut_str):
    match = re.match(r"([A-Z])(\d+)([A-Z])", mut_str.upper())
    if not match: return None
    wt, pos, mut = match.groups()
    pos = int(pos) - 1
    if pos >= len(seq): return None
    masked = seq[:pos] + "<mask>" + seq[pos+1:]
    inputs = tokenizer(masked, return_tensors="pt")
    with torch.no_grad():
        logits = model(**inputs).logits
    mask_pos = (inputs.input_ids == tokenizer.mask_token_id).nonzero()[0, 1]
    probs = torch.nn.functional.softmax(logits[0, mask_pos], dim=-1)
    wt_id = tokenizer.convert_tokens_to_ids(wt)
    mut_id = tokenizer.convert_tokens_to_ids(mut)
    return round(torch.log(probs[mut_id] / probs[wt_id]).item(), 3)

def get_pdb(uniprot_id):
    base = uniprot_id.split('-')[0]
    for v in ['v4', 'v3']:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{base}-F1-model_{v}.pdb"
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                return r.text, v
        except:
            continue
    return None, None

def render_3d(pdb_data, mutation):
    view = py3Dmol.view(width=1000, height=600)
    view.addModel(pdb_data, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.setBackgroundColor("#141e30")
    if mutation:
        pos = re.search(r"\d+", mutation).group()
        view.addStyle({"resi": int(pos)}, {"stick": {"color": "red", "radius": 0.5}})
        view.addStyle({"resi": int(pos)}, {"sphere": {"color": "red", "radius": 2, "opacity": 0.9}})
        view.zoomTo({"resi": int(pos)})
    view.zoom(1.3)
    view.spin(True)
    showmol(view, height=600, width=1000)

def get_plddt(pdb_data):
    scores = []
    for line in pdb_data.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                scores.append(float(line[60:66]))
            except:
                continue
    if scores:
        df = pd.DataFrame({"Residue": range(1, len(scores)+1), "pLDDT": scores})
        return df, round(sum(scores)/len(scores), 2)
    return None, None

# ========================= SIDEBAR =========================
with st.sidebar:
    st.markdown("# DNA GenoScope")
    st.markdown("**AI-Powered Mutation Intelligence**")
    st.divider()
    st.markdown("### Quick Start")
    st.code("UniProt: P04637\nMutation: R273H")
    st.divider()
    st.markdown("### Powered By")
    st.markdown("- DeepMind AlphaFold")
    st.markdown("- Meta AI ESM-2")
    st.markdown("- UniProt REST API")

# ========================= HEADER =========================
st.markdown("<h1 style='text-align:center; margin:2rem 0;'>GenoScope</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center; font-size:1.2rem; color:#aaa;'>Next-generation protein mutation analysis with AI & 3D visualization</p>", unsafe_allow_html=True)

# ========================= INPUT =========================
with st.container():
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown("### Enter Protein & Mutation")
        uniprot = st.text_input("UniProt ID", value="P04637", placeholder="e.g., P04637")
        mutation = st.text_input("Mutation", value="R273H", placeholder="e.g., R273H")
        
        col_a, col_b = st.columns(2)
        with col_a:
            analyze = st.button("Analyze Mutation", use_container_width=True, type="primary")
        with col_b:
            batch = st.file_uploader("Batch (.txt)", type="txt")

# ========================= ANALYSIS =========================
if analyze:
    with st.spinner("Analyzing with AI..."):
        bar = st.progress(0)
        status = st.empty()
        
        for i in range(100):
            time.sleep(0.015)
            bar.progress(i + 1)
        
        status.text("Fetching sequence...")
        seq, name = get_sequence(uniprot)
        
        if not seq:
            st.error("Invalid UniProt ID")
            st.stop()
            
        status.text("Computing mutation impact...")
        delta = compute_delta(seq, mutation)
        
        status.text("Loading 3D structure...")
        pdb, version = get_pdb(uniprot)
        
        status.text("Extracting confidence...")
        plddt_df, avg_plddt = get_plddt(pdb) if pdb else (None, None)
        
        impact = "Pathogenic" if delta and delta < -0.5 else "Likely Pathogenic" if delta and delta < 0 else "Benign"
        
        st.session_state.result = {
            "name": name,
            "mutation": mutation,
            "delta": delta,
            "impact": impact,
            "pdb": pdb,
            "version": version,
            "plddt_df": plddt_df,
            "avg_plddt": avg_plddt,
            "insight": INSIGHTS.get(mutation, "This mutation may significantly alter protein structure and function.")
        }
        
        celebrate()
        status.success("Analysis Complete!")

# ========================= RESULTS =========================
if "result" in st.session_state:
    r = st.session_state.result
    
    st.divider()
    st.markdown(f"### {r['name']} • {r['mutation']}")
    
    # Metrics Row
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.markdown(f"<div class='metric-card'><h3>Impact</h3><h2 style='color:#ff6b6b;'>{r['impact']}</h2>Δlog = {r['delta']}</div>", unsafe_allow_html=True)
    with c2:
        st.markdown(f"<div class='metric-card'><h3>pLDDT</h3><h2>{r['avg_plddt'] or 'N/A'}</h2></div>", unsafe_allow_html=True)
    with c3:
        st.markdown(f"<div class='metric-card'><h3>Model</h3><h2>{r['version']}</h2></div>", unsafe_allow_html=True)
    with c4:
        st.markdown(f"<div class='metric-card'><h3>Length</h3><h2>{len(seq)} aa</h2></div>", unsafe_allow_html=True)
    
    # 3D + Graph
    col_left, col_right = st.columns([65, 35])
    with col_left:
        st.markdown("### Interactive 3D Structure")
        if r['pdb']:
            render_3d(r['pdb'], r['mutation'])
            st.caption("Red sphere = Mutation site • Auto-rotating • Drag to interact")
    
    with col_right:
        st.markdown("### pLDDT Confidence Map")
        if r['plddt_df'] is not None:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.fill_between(r['plddt_df']['Residue'], r['plddt_df']['pLDDT'], color='#00c9ff', alpha=0.7)
            ax.plot(r['plddt_df']['Residue'], r['plddt_df']['pLDDT'], color='white', linewidth=2)
            ax.set_ylim(0, 100)
            ax.set_facecolor('#141e30')
            ax.spines['bottom'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.tick_params(colors='white')
            ax.yaxis.label.set_color('white')
            ax.xaxis.label.set_color('white')
            fig.patch.set_facecolor('#141e30')
            st.pyplot(fig)
        else:
            st.info("pLDDT data not available")
    
    # AI Insight
    st.markdown("### AI Clinical Insight")
    st.markdown(f"<div class='insight-box'>{r['insight']}</div>", unsafe_allow_html=True)

else:
    st.markdown("""
    <div style='text-align:center; padding:4rem; color:#888;'>
        <h3>Ready to explore cancer-causing mutations?</h3>
        <p>Use P04637 + R273H to see full analysis with 3D spinning protein</p>
    </div>
    """, unsafe_allow_html=True)

# ========================= FOOTER =========================
st.markdown("""
<div style='text-align:center; padding:2rem; color:#666; margin-top:4rem;'>
    © 2025 GenoScope • Built with Streamlit • Powered by AlphaFold & ESM-2
</div>
""", unsafe_allow_html=True)