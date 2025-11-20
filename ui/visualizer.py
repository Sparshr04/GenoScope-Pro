# ui/visualizer.py
import py3Dmol
from stmol import showmol
import matplotlib.pyplot as plt
import pandas as pd
import re
import streamlit as st

def render_3d_structure(pdb_data: str, mutation: str = None):
    """Renders the interactive 3D molecule."""
    if not pdb_data:
        return

    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()

    # Highlight mutation
    if mutation:
        m = re.match(r"[A-Z](\d+)[A-Z]", mutation.upper())
        if m:
            pos = int(m.group(1))
            # Highlight specific residue
            view.addStyle({"resi": str(pos)}, {"stick": {"color": "red", "radius": 0.3}})
            view.zoomTo({"resi": str(pos)})
    
    showmol(view, height=600, width=800)

def plot_plddt_confidence(pdb_data: str):
    """Extracts B-factors (pLDDT) and plots confidence graph."""
    if not pdb_data:
        return

    # Parsing logic (Simplified for brevity)
    lines = pdb_data.splitlines()
    res_plddt = {}
    for line in lines:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            resnum = int(line[22:26].strip())
            plddt = float(line[60:66])
            res_plddt[resnum] = plddt
            
    if not res_plddt:
        return

    df = pd.DataFrame(list(res_plddt.items()), columns=['Residue', 'pLDDT'])
    
    # Matplotlib Plot
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(df["Residue"], df["pLDDT"], color='blue', linewidth=1)
    ax.axhline(y=90, color='green', linestyle='--', alpha=0.5)
    ax.axhline(y=70, color='orange', linestyle='--', alpha=0.5)
    ax.axhline(y=50, color='red', linestyle='--', alpha=0.5)
    ax.set_title("AlphaFold Confidence (pLDDT) per Residue")
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("Confidence Score")
    
    st.pyplot(fig)
