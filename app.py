import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# 1. Page Configuration
st.set_page_config(page_title="Bioinformatics 3D Portal", page_icon="🧬", layout="wide")

st.title("🧬 Next-Gen 3D Protein Structure & Analytics Platform")
st.markdown("Visualize 3D protein geometry and calculate primary sequence stability metrics.")
st.markdown("---")

# 2. Sidebar Inputs
st.sidebar.header("📁 Step 1: Load 3D Structure")
input_method = st.sidebar.radio(
    "Choose Structure Source:",
    ["Use Built-in Demo Protein (GFP - 1EMA)", "Upload .pdb File from Colab"]
)

pdb_string = None

if input_method == "Upload .pdb File from Colab":
    uploaded_file = st.sidebar.file_uploader("Upload .pdb file", type=["pdb", "txt"])
    if uploaded_file is not None:
        raw_text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
        if "ATOM" in raw_text or "HETATM" in raw_text:
            pdb_string = raw_text
        else:
            st.sidebar.error("⚠️ Uploaded file has no valid 3D ATOM coordinates. Loading backup GFP structure automatically!")
            pdb_string = requests.get("https://rcsb.org").text
else:
    with st.spinner("Fetching high-resolution 3D coordinates (1EMA)..."):
        pdb_string = requests.get("https://rcsb.org").text

st.sidebar.header("📝 Step 2: Sequence Analytics")
default_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
sequence_input = st.sidebar.text_area("Paste Raw Amino Acid Sequence:", default_seq, height=150)

# 3. Main Dashboard Logic
if pdb_string:
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("👁️ Interactive 3D Molecular Canvas")
        
        c1, c2 = st.columns(2)
        render_style = c1.selectbox("Render Style:", ["cartoon", "stick", "sphere"])
        color_scheme = c2.selectbox("Color Theme:", ["spectrum", "chain", "b-factor"])
        
        # Build 3D Viewer using native py3Dmol HTML generator
        view = py3Dmol.view(width=600, height=450)
        view.addModel(pdb_string, 'pdb')
        
        if color_scheme == "b-factor":
            view.setStyle({render_style: {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}})
        else:
            view.setStyle({render_style: {'color': color_scheme}})
            
        view.setBackgroundColor('white')
        view.zoomTo()
        
        # Render directly inside Streamlit iframe
        components.html(view._make_html(), height=470, width=620)
        st.caption("💡 Tip: Left-click + drag to rotate, scroll to zoom, right-click + drag to pan.")
        
    with col2:
        st.subheader("📊 Primary Sequence & Stability Metrics")
        
        if sequence_input:
            clean_seq = "".join([line.strip() for line in sequence_input.splitlines() if not line.startswith(">")])
            clean_seq = "".join([c for c in clean_seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY"])
            
            if len(clean_seq) > 0:
                analysed_seq = ProteinAnalysis(clean_seq)
                mw = analysed_seq.molecular_weight()
                aromaticity = analysed_seq.aromaticity()
                instability = analysed_seq.instability_index()
                isoelectric = analysed_seq.isoelectric_point()
                
                m1, m2 = st.columns(2)
                m1.metric("Molecular Weight", f"{mw:,.2f} Da")
                m2.metric("Isoelectric Point (pI)", f"{isoelectric:.2f}")
                
                m3, m4 = st.columns(2)
                m3.metric("Aromaticity", f"{aromaticity:.3f}")
                m4.metric("Instability Index", f"{instability:.2f}")
                
                st.markdown("#### Structural Stability Prediction")
                if instability < 40:
                    st.success("✅ **STABLE:** Instability index is below 40.0. The protein is predicted to be stable in a test tube.")
                else:
                    st.warning("⚠️ **UNSTABLE:** Instability index is above 40.0. The protein may degrade quickly without chaperones.")
            else:
                st.error("⚠️ Please enter valid standard amino acid letters.")
        
        st.markdown("---")
        st.subheader("📥 Data Export Vault")
        st.download_button(
            label="💾 Download Verified .PDB File",
            data=pdb_string,
            file_name="verified_structure.pdb",
            mime="text/plain"
        )
else:
    st.info("👈 **Waiting for Data:** Select the Built-in Demo Protein or upload a `.pdb` file in the left sidebar!")
