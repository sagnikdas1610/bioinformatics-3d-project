import streamlit as st
import streamlit.components.v1 as components
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json

# 1. Page Configuration
st.set_page_config(page_title="Bioinformatics 3D Portal", page_icon="🧬", layout="wide")

st.title("🧬 Next-Gen 3D Protein Structure & Analytics Platform")
st.markdown("Upload a `.pdb` file generated from our Google Colab AI Engine to visualize its 3D geometry and calculate stability metrics.")
st.markdown("---")

# 2. Sidebar Inputs
st.sidebar.header("📁 Step 1: Upload 3D Structure")
uploaded_file = st.sidebar.file_uploader("Upload .pdb file from Colab", type=["pdb"])

st.sidebar.header("📝 Step 2: Sequence Analytics")
default_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
sequence_input = st.sidebar.text_area("Paste Raw Amino Acid Sequence (No FASTA header):", default_seq, height=150)

# 3. Direct 3Dmol HTML Renderer (Bypasses stmol loading errors)
def render_3d_molecule(pdb_data, style_type, color_theme):
    pdb_json = json.dumps(pdb_data)
    
    if color_theme == "b-factor":
        style_str = f"{{{style_type}: {{colorscheme: {{prop: 'b', gradient: 'roygb', min: 50, max: 90}}}}}}"
    else:
        style_str = f"{{{style_type}: {{color: '{color_theme}'}}}}"

    html_code = f"""
    <script src="https://cloudflare.com"></script>
    <div id="container-01" style="width: 100%; height: 450px; position: relative; border-radius: 8px; overflow: hidden;"></div>
    <script>
        let element = document.querySelector('#container-01');
        let config = {{ backgroundColor: '#ffffff' }};
        let viewer = $3Dmol.createViewer(element, config);
        viewer.addModel({pdb_json}, "pdb");
        viewer.setStyle({{}}, {style_str});
        viewer.zoomTo();
        viewer.render();
    </script>
    """
    components.html(html_code, height=460)

# 4. Main Dashboard Logic
if uploaded_file is not None:
    pdb_string = uploaded_file.getvalue().decode("utf-8")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("👁️ Interactive 3D Molecular Canvas")
        
        c1, c2 = st.columns(2)
        render_style = c1.selectbox("Render Style:", ["cartoon", "stick", "sphere"])
        color_scheme = c2.selectbox("Color Theme:", ["spectrum", "chain", "b-factor"])
        
        # Render the 3D Molecule directly
        render_3d_molecule(pdb_string, render_style, color_scheme)
        st.caption("💡 Tip: Use your mouse to rotate (left-click), zoom (scroll), and pan (right-click).")
        
    with col2:
        st.subheader("📊 Primary Sequence & Stability Metrics")
        
        if sequence_input:
            clean_seq = "".join([line.strip() for line in sequence_input.splitlines() if not line.startswith(">")])
            clean_seq = clean_seq.replace(" ", "").upper()
            
            try:
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
                    
            except Exception as e:
                st.error("⚠️ Invalid characters detected in sequence. Ensure only standard single-letter amino acid codes are used.")
        
        st.markdown("---")
        st.subheader("📥 Data Export Vault")
        st.download_button(
            label="💾 Download Verified .PDB File",
            data=pdb_string,
            file_name="verified_structure.pdb",
            mime="text/plain"
        )
else:
    st.info("👈 **Waiting for Data:** Please run the Google Colab Engine to generate your `.pdb` structure, then upload it using the sidebar on the left!")
