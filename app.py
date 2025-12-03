import streamlit as st
import requests
import py3Dmol
import pandas as pd
import io

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED
from Bio.PDB import PDBParser

# ----------------------------------------------------
# Page setup
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# Minimal CSS (Just for scrollbars & headers)
# ----------------------------------------------------
st.markdown(
    """
    <style>
    /* Remove top padding to make header flush */
    .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }

    /* Style for the metadata scroll area */
    .meta-scroll {
        max-height: 75vh;
        overflow-y: auto;
        padding-right: 10px;
    }

    /* Scrollbar styling for Webkit */
    .meta-scroll::-webkit-scrollbar {
        width: 6px;
    }
    .meta-scroll::-webkit-scrollbar-track {
        background: transparent;
    }
    .meta-scroll::-webkit-scrollbar-thumb {
        background-color: #ccc;
        border-radius: 20px;
    }

    /* Subtle divider */
    .sub-divider {
        margin-top: 1rem;
        margin-bottom: 1rem;
        border-top: 1px solid #f0f2f6;
    }
    
    /* Text styling */
    .label-text {
        font-weight: 600;
        font-size: 0.8rem;
        color: #666;
        margin-bottom: 0px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    .value-text {
        font-size: 0.95rem;
        color: #111;
        margin-bottom: 12px;
        word-wrap: break-word;
        font-family: 'Source Code Pro', monospace;
    }
    
    .metric-container {
        background-color: #f8f9fa;
        padding: 10px;
        border-radius: 8px;
        margin-bottom: 10px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# Constants & URLs
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents/PDBs"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/PDBs/"
METADATA_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_data_2811.xlsx"
LOGO_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs.png"

# ----------------------------------------------------
# Data Fetching Functions
# ----------------------------------------------------
@st.cache_data
def list_pdb_files():
    try:
        r = requests.get(GITHUB_API_URL)
        r.raise_for_status()
        files = r.json()
        pdb_files = [
            f["name"]
            for f in files
            if isinstance(f, dict) and f.get("name", "").lower().endswith(".pdb")
        ]
        return sorted(pdb_files)
    except Exception as e:
        st.error(f"Error fetching file list: {e}")
        return []

def fetch_pdb_from_github(filename: str) -> str | None:
    try:
        url = f"{GITHUB_RAW_BASE}{filename}"
        r = requests.get(url)
        if r.status_code == 200 and r.text.strip():
            return r.text
        return None
    except:
        return None

@st.cache_data
def load_metadata():
    try:
        df = pd.read_excel(METADATA_URL)
        df.columns = [c.strip().lower() for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

# ----------------------------------------------------
# Computation Functions
# ----------------------------------------------------
def calculate_esol(mol, logp, mw, rb, aromatic_rings):
    """
    Estimate solubility (ESOL)
    LogS = 0.16 - 0.63(cLogP) - 0.0062(MW) + 0.066(RB) - 0.74(AromaticProportion)
    """
    try:
        num_heavy = mol.GetNumHeavyAtoms()
        num_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_prop = num_aromatic / num_heavy if num_heavy > 0 else 0
        
        # ESOL Formula
        esol = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * aromatic_prop)
        return esol
    except:
        return None

def compute_physchem_from_pdb(pdb_text: str) -> dict:
    props = {}
    try:
        mol = Chem.MolFromPDBBlock(pdb_text, sanitize=True, removeHs=False)
        if mol is None:
            return props

        # Basic Descriptors
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rb = Lipinski.NumRotatableBonds(mol)
        
        # Advanced Descriptors
        qed = QED.qed(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        
        # ESOL Calculation
        esol = calculate_esol(mol, logp, mw, rb, aromatic_rings)

        # Store nicely formatted strings
        props["Mol Wt"] = f"{mw:.2f}"
        props["LogP"] = f"{logp:.2f}"
        props["TPSA"] = f"{rdMolDescriptors.CalcTPSA(mol):.2f}"
        props["QED"] = f"{qed:.3f}"
        props["ESOL (LogS)"] = f"{esol:.2f}" if esol else "N/A"
        
        props["H-Acc"] = Lipinski.NumHAcceptors(mol)
        props["H-Don"] = Lipinski.NumHDonors(mol)
        props["Rot. Bonds"] = rb
        props["Arom. Rings"] = aromatic_rings
        props["Sat. Rings"] = Lipinski.NumSaturatedRings(mol)
        props["Atoms"] = mol.GetNumAtoms()
        props["F-Csp3"] = f"{rdMolDescriptors.CalcFractionCSP3(mol):.2f}"
        
    except Exception:
        pass
    return props

def compute_biopython_features(pdb_text: str) -> dict:
    feats = {}
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("lig", io.StringIO(pdb_text))
        # Just simple counts
        feats["Bio Res"] = sum(1 for _ in structure.get_residues())
        feats["Bio Atoms"] = sum(1 for _ in structure.get_atoms())
    except Exception:
        pass
    return feats

def find_metadata(metadata_df, pdb_filename):
    if metadata_df.empty:
        return None
        
    pdb_root = pdb_filename.replace(".pdb", "").lower()
    if "pdbs" not in metadata_df.columns:
        return None
    
    metadata_df["match"] = metadata_df["pdbs"].astype(str).str.lower()
    
    # Try exact match
    match = metadata_df[metadata_df["match"] == pdb_filename.lower()]
    if not match.empty:
        return match

    # Try root match
    match = metadata_df[metadata_df["match"] == pdb_root]
    return match if not match.empty else None

# ----------------------------------------------------
# 3D Viewer Function
# ----------------------------------------------------
def show_3d_pdb(pdb_text: str):
    # Increased height for better immersion
    view = py3Dmol.view(width="100%", height=700)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    view.setBackgroundColor("white") # Matches streamlit light theme usually
    html = view._make_html()
    st.components.v1.html(html, height=700)

# ----------------------------------------------------
# APP LOGIC
# ----------------------------------------------------

# 1. SIDEBAR CONTROLS
with st.sidebar:
    st.image(LOGO_URL, use_container_width=True)
    st.markdown("### NucLigs Database")
    
    # Load Data
    all_pdb_files = list_pdb_files()
    metadata_df = load_metadata()

    # Search & Select
    search_query = st.text_input("🔍 Search Structure", placeholder="e.g., 1A2B...")
    
    if search_query:
        pdb_files = [p for p in all_pdb_files if search_query.lower() in p.lower()]
        if not pdb_files:
            st.warning("No matches found.")
            pdb_files = all_pdb_files
    else:
        pdb_files = all_pdb_files

    selected_pdb = st.selectbox("Select PDB File", pdb_files, index=0 if pdb_files else None)
    
    st.markdown("---")
    st.caption(f"**Source:** tushar1298/qwertyui")
    st.caption(f"**Total Files:** {len(all_pdb_files)}")
    if selected_pdb:
         st.caption(f"**Viewing:** {selected_pdb}")

# 2. MAIN AREA
if not selected_pdb:
    st.info("Please select a structure from the sidebar to begin.")
else:
    pdb_text = fetch_pdb_from_github(selected_pdb)

    if pdb_text:
        # Compute properties
        physchem = compute_physchem_from_pdb(pdb_text)
        
        # Main Layout: 3 Columns
        # 1. Viewer (Largest)
        # 2. Predicted Features (Medium)
        # 3. Metadata (Medium)
        col_viewer, col_preds, col_meta = st.columns([2, 0.8, 0.8])

        with col_viewer:
            st.subheader(f"Structure Visualization")
            show_3d_pdb(pdb_text)

        # --- NEW COLUMN: PREDICTED FEATURES ---
        with col_preds:
            st.subheader("Predicted Features")
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            
            if physchem:
                st.markdown("##### 🧪 Physical")
                st.markdown(f"""
                <div class="metric-container">
                    <div class="label-text">Est. Solubility (ESOL)</div>
                    <div class="value-text">{physchem.get("ESOL (LogS)", "-")}</div>
                    <div class="label-text">LogP</div>
                    <div class="value-text">{physchem.get("LogP", "-")}</div>
                    <div class="label-text">Mol. Weight</div>
                    <div class="value-text">{physchem.get("Mol Wt", "-")} g/mol</div>
                </div>
                """, unsafe_allow_html=True)

                st.markdown("##### 💊 Druglikeness")
                st.markdown(f"""
                <div class="metric-container">
                    <div class="label-text">QED Score</div>
                    <div class="value-text">{physchem.get("QED", "-")}</div>
                    <div class="label-text">TPSA</div>
                    <div class="value-text">{physchem.get("TPSA", "-")} Å²</div>
                    <div class="label-text">Fraction Csp3</div>
                    <div class="value-text">{physchem.get("F-Csp3", "-")}</div>
                </div>
                """, unsafe_allow_html=True)

                st.markdown("##### 🧬 Composition")
                c1, c2 = st.columns(2)
                c1.metric("H-Acc", physchem.get("H-Acc", "-"))
                c2.metric("H-Don", physchem.get("H-Don", "-"))
                c1.metric("Rot Bonds", physchem.get("Rot. Bonds", "-"))
                c2.metric("Atoms", physchem.get("Atoms", "-"))
                c1.metric("Arom Rings", physchem.get("Arom. Rings", "-"))
                c2.metric("Sat Rings", physchem.get("Sat. Rings", "-"))
            
            else:
                st.warning("Could not compute properties.")
            st.markdown("</div>", unsafe_allow_html=True)


        # --- COLUMN: REPOSITORY METADATA ---
        with col_meta:
            st.subheader("Metadata")
            st.markdown('<div class="meta-scroll">', unsafe_allow_html=True)
            
            row = find_metadata(metadata_df, selected_pdb)

            if row is not None and not row.empty:
                data = row.iloc[0].to_dict()
                data.pop("match", None)
                
                long_fields = ["names", "smiles", "inchi", "description"]
                
                # Standard fields
                st.markdown("##### 📋 Details")
                for key, value in data.items():
                    if key.lower() not in long_fields:
                        st.markdown(
                            f"""
                            <div class="label-text">{key.replace('_', ' ').title()}</div>
                            <div class="value-text">{value}</div>
                            """, 
                            unsafe_allow_html=True
                        )
                
                st.markdown('<div class="sub-divider"></div>', unsafe_allow_html=True)

                # Long fields
                for key in long_fields:
                    if key in data:
                        st.markdown(f'<div class="label-text">{key.upper()}</div>', unsafe_allow_html=True)
                        st.code(str(data[key]), language="text")
            else:
                st.info("No external metadata found.")
            
            st.markdown("</div>", unsafe_allow_html=True)
