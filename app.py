# app.py
import streamlit as st
from core.ai_engine import MutationPredictor
from core.data_client import BioDataClient
from ui.visualizer import render_3d_structure, plot_plddt_confidence

st.set_page_config(page_title="GenoScope Pro", page_icon="🧬", layout="wide")


# These are lightweight classes now; heavy loading happens inside via caching
predictor = MutationPredictor()
client = BioDataClient()

# --- UI LAYOUT ---
st.title("🧬 GenoScope Pro")
st.markdown("Modular AI Mutation Analysis System")
st.divider()

col1, col2 = st.columns([2, 1])

with col1:
    uniprot_id = st.text_input("UniProt ID", value="P04637")
    mutation = st.text_input("Mutation", value="R273H")
    gene_name = st.text_input("Gene Name", value="TP53")
    
    if st.button("🚀 Run Analysis"):
        with st.spinner("Processing..."):
            # 1. Fetch Data
            seq, desc = client.get_sequence(uniprot_id)
            
            if not seq:
                st.error("Could not fetch sequence.")
            else:
                # 2. Run AI
                try:
                    score = predictor.compute_score(seq, mutation)
                    result = predictor.analyze_impact(score, mutation, gene_name)
                    
                    # Store in session for other tabs
                    st.session_state['result'] = result
                    st.session_state['seq'] = seq
                    st.session_state['desc'] = desc
                except ValueError as e:
                    st.error(str(e))

# --- RESULTS TABS ---
tab1, tab2, tab3 = st.tabs(["🔍 AI Analysis", "🧬 3D Structure", "🏥 Clinical Data"])

with tab1:
    if 'result' in st.session_state:
        res = st.session_state['result']
        st.subheader(f"Verdict: {res['severity']}")
        st.info(res['insight'])
        st.caption(f"Source: {res['source']} | Raw Delta: {res['score']:.4f}")
        st.text_area("Sequence", st.session_state['seq'], height=100)

with tab2:
    if st.button("Load 3D Model"):
        pdb_text, ver = client.get_structure(uniprot_id)
        if pdb_text:
            st.success(f"Loaded AlphaFold Model ({ver})")
            col_a, col_b = st.columns([3, 1])
            with col_a:
                render_3d_structure(pdb_text, mutation)
            with col_b:
                plot_plddt_confidence(pdb_text)
        else:
            st.warning("Structure not found.")

with tab3:
    if st.button("🏥 Check Clinical Databases"):
        with st.spinner(f"Searching ClinVar for {gene_name} {mutation}..."):
            clinical_data = client.fetch_clinical_data(gene_name, mutation)
            
            if "error" in clinical_data:
                st.error(f"Error: {clinical_data['error']}")
            elif clinical_data.get("significance") == "Not Found":
                st.warning("No direct match found in ClinVar for this specific mutation.")
                st.caption("Try checking the gene name or mutation format.")
            else:
                st.success("Record Found!")
                col_a, col_b = st.columns(2)
                with col_a:
                    st.metric("Clinical Significance", clinical_data.get("significance"))
                with col_b:
                    st.write("**Associated Conditions:**")
                    st.write(clinical_data.get("conditions"))
                st.caption(f"Source: {clinical_data.get('source')}")

