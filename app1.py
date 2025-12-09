import streamlit as st
import requests
import py3Dmol
import pandas as pd
import io
import zipfile

from supabase import create_client, Client

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED
from Bio.PDB import PDBParser

# ----------------------------------------------------
# Page setup (Must be first)
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# CSS Styling
# ----------------------------------------------------
st.markdown(
    """
    <style>
    /* Global Container Padding */
    .block-container {
        padding-top: 3.5rem;
        padding-bottom: 3rem;
    }

    /* Scrollable Containers */
    .meta-scroll {
        max-height: 70vh;
        overflow-y: auto;
        padding-right: 12px;
    }
    .meta-scroll::-webkit-scrollbar {
        width: 8px;
    }
    .meta-scroll::-webkit-scrollbar-track {
        background: #f1f1f1;
        border-radius: 4px;
    }
    .meta-scroll::-webkit-scrollbar-thumb {
        background: #c1c1c1;
        border-radius: 4px;
    }
    .meta-scroll::-webkit-scrollbar-thumb:hover {
        background: #a8a8a8;
    }

    /* Section Cards in DB View */
    .feature-card {
        background-color: #ffffff;
        border: 1px solid #e0e0e0;
        border-radius: 10px;
        padding: 15px;
        margin-bottom: 15px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.03);
    }
    .feature-card h5 {
        color: #2c3e50;
        font-size: 0.95rem;
        font-weight: 700;
        margin-bottom: 12px;
        border-bottom: 2px solid #f0f2f6;
        padding-bottom: 8px;
    }
    
    /* Reference Card specific styling */
    .ref-card {
        background-color: #fcfcfc;
        border-left: 4px solid #3498db;
        padding: 15px;
        margin-bottom: 15px;
        border-radius: 6px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .ref-title {
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 8px;
        font-size: 1.0rem;
        line-height: 1.4;
    }
    .ref-meta {
        font-size: 0.85rem;
        color: #555;
        margin-bottom: 8px;
    }

    /* ID Highlight Card */
    .id-card {
        background-color: #f8f9fa;
        border-left: 5px solid #4CAF50;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    .id-label {
        font-size: 0.75rem;
        color: #666;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .id-value {
        font-size: 1.5rem;
        font-weight: 800;
        color: #2c3e50;
        margin: 4px 0 8px 0;
        line-height: 1.2;
    }
    .id-sub {
        font-size: 1.0rem;
        color: #555;
        font-weight: 500;
        font-style: italic;
    }

    /* Data Points */
    .data-row {
        display: flex;
        justify-content: space-between;
        margin-bottom: 8px;
        align-items: flex-start;
    }

    /* Tighter layout for rows inside reference cards */
    .ref-card .data-row {
        justify-content: flex-start;
        align-items: center;
    }
    .ref-card .data-label {
        min-width: 55px;
        margin-right: 2px;
    }
    .ref-card .data-value {
        text-align: left;
    }

    .data-label {
        font-size: 0.85rem;
        color: #666;
        font-weight: 600;
        min-width: 140px;
        flex-shrink: 0;
        margin-right: 15px;
        padding-top: 2px;
    }
    .data-value {
        font-family: 'Source Code Pro', monospace;
        font-size: 0.9rem;
        color: #222;
        font-weight: 600;
        text-align: right;
        word-break: break-word;
        flex-grow: 1;
    }
    .reference-text {
        font-size: 0.75rem;
        color: #999;
        margin-left: 5px;
        font-weight: 400;
    }

    /* Homepage Cards */
    .home-card {
        padding: 25px;
        border-radius: 12px;
        border: 1px solid #eef2f5;
        background-color: white;
        text-align: left;
        transition: all 0.2s ease;
        box-shadow: 0 4px 6px rgba(0,0,0,0.02);
        height: 100%;
    }
    .home-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 12px 20px rgba(0,0,0,0.08);
        border-color: #4CAF50;
    }
    .home-card h3 {
        color: #2c3e50;
        font-size: 1.2rem;
        margin-bottom: 15px;
        font-weight: 700;
    }
    .home-card p {
        color: #555;
        font-size: 0.95rem;
        line-height: 1.6;
    }
    
    /* Buttons */
    .stButton > button {
        width: 100%;
        border-radius: 8px;
        font-weight: 600;
        padding-top: 0.5rem;
        padding-bottom: 0.5rem;
    }
    
    /* Tabs Styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 10px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: #f8f9fa;
        border-radius: 8px;
        color: #495057;
        font-weight: 600;
        padding: 0 20px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #e8f5e9;
        color: #2e7d32;
        border-bottom: 2px solid #2e7d32;
    }
    
    /* Badge Styles */
    .badge-pass { background-color: #d4edda; color: #155724; padding: 4px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 700; }
    .badge-fail { background-color: #f8d7da; color: #721c24; padding: 4px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 700; }
    
    /* Sidebar Logo Adjustment */
    [data-testid="stSidebar"] img {
        margin-bottom: 0px; 
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------
# Supabase Configuration
# ----------------------------------------------------
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"

BUCKET_NAME = "NucLigs_PDBs"       # PDB files
METADATA_BUCKET = "codes"          # Excel files
METADATA_FILENAME = "NucLigs_metadata.xlsx"
METADATA_REF_FILENAME = "references.xlsx"

@st.cache_resource
def init_supabase():
    try:
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception:
        return None

supabase = init_supabase()

# ----------------------------------------------------
# External URLs
# ----------------------------------------------------
LOGO_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs.png"

# ----------------------------------------------------
# Data & Compute Functions
# ----------------------------------------------------
@st.cache_data(ttl=0)
def load_metadata():
    """Load main metadata sheet, normalize column names."""
    if not supabase:
        return pd.DataFrame()
    try:
        data_bytes = supabase.storage.from_(METADATA_BUCKET).download(METADATA_FILENAME)
        df = pd.read_excel(io.BytesIO(data_bytes))
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
        return df
    except Exception as e:
        st.sidebar.error(f"Error loading metadata from Supabase: {e}")
        return pd.DataFrame()

@st.cache_data(ttl=0)
def load_references():
    """Load reference sheet, normalize columns."""
    if not supabase:
        return pd.DataFrame()
    try:
        data_bytes = supabase.storage.from_(METADATA_BUCKET).download(METADATA_REF_FILENAME)
        df = pd.read_excel(io.BytesIO(data_bytes), sheet_name=0)
        df.columns = [str(c).strip().lower().replace(" ", "_") for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()

@st.cache_data(ttl=0)
def get_ids_from_metadata():
    df = load_metadata()
    if not df.empty and 'nl' in df.columns:
        return sorted(df['nl'].dropna().astype(str).unique().tolist())
    return []

def fetch_pdb_from_supabase(filename_or_id: str) -> str | None:
    if not supabase:
        return None
    try:
        filename = filename_or_id if filename_or_id.lower().endswith(".pdb") else f"{filename_or_id}.pdb"
        data_bytes = supabase.storage.from_(BUCKET_NAME).download(filename)
        return data_bytes.decode('utf-8')
    except Exception:
        return None

def calculate_esol(mol, logp, mw, rb, aromatic_rings):
    try:
        num_heavy = mol.GetNumHeavyAtoms()
        num_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_prop = num_aromatic / num_heavy if num_heavy > 0 else 0
        esol = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * aromatic_prop)
        return esol
    except Exception:
        return None

def compute_physchem(mol) -> dict:
    props = {}
    if mol is None:
        return props
    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_acc = Lipinski.NumHAcceptors(mol)
        h_don = Lipinski.NumHDonors(mol)
        rb = Lipinski.NumRotatableBonds(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        charge = Chem.GetFormalCharge(mol)
        chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        qed = QED.qed(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        esol = calculate_esol(mol, logp, mw, rb, aromatic_rings)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if h_don > 5: violations += 1
        if h_acc > 10: violations += 1

        props["Formula"] = formula
        props["Charge"] = str(charge)
        props["Chiral Centers"] = str(chiral_centers)
        props["Mol Wt"] = f"{mw:.2f}"
        props["LogP"] = f"{logp:.2f}"
        props["TPSA"] = f"{rdMolDescriptors.CalcTPSA(mol):.2f}"
        props["QED"] = f"{qed:.3f}"
        props["ESOL (LogS)"] = f"{esol:.2f}" if esol is not None else "N/A"
        props["H-Acc"] = h_acc
        props["H-Don"] = h_don
        props["Rot. Bonds"] = rb
        props["Arom. Rings"] = aromatic_rings
        props["Sat. Rings"] = Lipinski.NumSaturatedRings(mol)
        props["Atoms"] = mol.GetNumAtoms()
        props["F-Csp3"] = f"{rdMolDescriptors.CalcFractionCSP3(mol):.2f}"
        props["Lipinski Violations"] = violations
        props["_RDKitMol"] = mol
    except Exception:
        pass
    return props

def show_3d_pdb(pdb_text: str, style_choice: str = "Stick", bg_color: str = "white"):
    """Render viewer and return the py3Dmol view object (for PNG export)."""
    view = py3Dmol.view(width=900, height=700)
    view.addModel(pdb_text, "pdb")
    
    if style_choice == "Stick":
        view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    elif style_choice == "Sphere":
        view.setStyle({"sphere": {"colorscheme": "greenCarbon"}})
    elif style_choice == "Line":
        view.setStyle({"line": {"colorscheme": "greenCarbon"}})
    else:
        view.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    
    view.zoomTo()
    view.setBackgroundColor(bg_color)
    html = view._make_html()
    st.components.v1.html(html, height=700)
    return view

def render_row(label, value, ref=None, help_text=None):
    ref_html = f'<span class="reference-text">({ref})</span>' if ref else ''
    st.markdown(
        f"""<div class="data-row">
               <span class="data-label">{label}</span>
               <span class="data-value">{value}{ref_html}</span>
           </div>""",
        unsafe_allow_html=True
    )

# ----------------------------------------------------
# PAGE RENDERERS
# ----------------------------------------------------
def render_homepage():
    st.markdown(
        f"""
        <div style="text-align: center; padding-top: 10px;">
            <img src="{LOGO_URL}" width="180" style="margin-bottom: 15px;">
            <h1 style='color: #2c3e50; margin-bottom: 0;'>NucLigs Database</h1>
            <p style='color: #666; font-size: 1.15rem; font-weight: 300;'>
                The Premier Resource for Nucleoside and Nucleotide analog Structures
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("""
    <div style='background-color: #f8f9fa; padding: 30px; border-radius: 12px; border-left: 5px solid #4CAF50; box-shadow: 0 4px 15px rgba(0,0,0,0.05); margin-bottom: 40px;'>
        <h3 style='color: #2c3e50; margin-top: 0;'>About the Database</h3>
        <p style='color: #444; font-size: 1.05rem; line-height: 1.6;'>
            The <b>NucLigs Database</b> is a specialized repository designed to facilitate research in the field of nucleic acid targeting by incorporating 
            <b>nucleoside and nucleotide analogs</b> at one place providing a unified platform for detailed analysis and visualization.
        </p>
        <p style='color: #444; font-size: 1.05rem; line-height: 1.6; margin-bottom: 0;'>
            Whether you are involved in <b>rational drug design</b>, <b>structural biology</b>, or <b>cheminformatics</b>, NucLigs offers robust tools to explore 
            the physico-chemical landscape of nucleic acid interactions, supporting the discovery of next-generation therapeutics.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    f1, f2, f3 = st.columns(3)
    with f1:
        st.markdown("""
        <div class="home-card">
            <h3>3D Visualization</h3>
            <p>Interactive, high-fidelity rendering of ligand-target complexes using Py3Dmol. Inspect binding modes, molecular surfaces, and structural conformations in real-time directly within your browser.</p>
        </div>
        """, unsafe_allow_html=True)
    with f2:
        st.markdown("""
        <div class="home-card">
            <h3>Chemical Profiling</h3>
            <p>Automated calculation of critical molecular descriptors. Access data on Molecular Weight, LogP, TPSA, and Lipinski's Rule of 5 compliance powered by the RDKit cheminformatics engine.</p>
        </div>
        """, unsafe_allow_html=True)
    with f3:
        st.markdown("""
        <div class="home-card">
            <h3>Data Accessibility</h3>
            <p>Seamlessly retrieve standardized structural data. Export ligands and complexes in industry-standard formats (PDB, SDF, MOL2) to integrate directly with your local modeling workflows.</p>
        </div>
        """, unsafe_allow_html=True)
        
    st.markdown("<br><br>", unsafe_allow_html=True)

    _, btn_col, _ = st.columns([1.5, 1, 1.5])
    with btn_col:
        if st.button("Explore the Collection", type="primary", use_container_width=True):
            st.session_state['page'] = 'database'
            st.rerun()
    
    st.markdown(
        "<div style='text-align: center; margin-top: 50px; color: #aaa; font-size: 0.85rem;'>© 2024 NucLigs Database Project • Version 2.2</div>",
        unsafe_allow_html=True
    )

def render_database():
    metadata_df = load_metadata()
    refs_df = load_references()
    all_nuc_ids = get_ids_from_metadata()
    
    # Map ID -> "Name (ID)"
    id_map = {}
    if not metadata_df.empty:
        temp_df = metadata_df.copy()
        if 'nl' in temp_df.columns:
            temp_df.set_index('nl', inplace=True)
            name_col = 'names' if 'names' in temp_df.columns else ('name' if 'name' in temp_df.columns else None)
            if name_col:
                for nid, row in temp_df.iterrows():
                    chem_name = str(row[name_col]) if pd.notna(row[name_col]) else "Unknown"
                    if len(chem_name) > 50:
                        chem_name = chem_name[:47] + "..."
                    id_map[str(nid)] = f"{chem_name} ({nid})"

    # Sidebar
    with st.sidebar:
        if st.button("Back to Home"):
            st.session_state['page'] = 'home'
            st.rerun()
            
        st.markdown("---")
        st.markdown("### Structure Finder")

        search_query = st.text_input(
            "Filter database:",
            placeholder="Search NucL ID or Name...",
            label_visibility="collapsed"
        )
        
        is_searching = False
        if search_query:
            is_searching = True
            q = search_query.lower()
            mask = (
                metadata_df['nl'].astype(str).str.lower().str.contains(q, na=False) |
                metadata_df['names'].astype(str).str.lower().str.contains(q, na=False)
            )
            filtered_matches = metadata_df[mask]['nl'].dropna().unique().tolist()
            nuc_ids = sorted([str(x) for x in filtered_matches])
            
            if not nuc_ids:
                st.warning("No matches found.")
                nuc_ids = []
            else:
                st.success(f"Found {len(nuc_ids)} structures")
        else:
            nuc_ids = all_nuc_ids

        if nuc_ids:
            format_strategy = (lambda x: id_map.get(x, x)) if is_searching else (lambda x: x)
            selected_nuc_id = st.selectbox(
                "Select Structure Result:",
                nuc_ids,
                index=0,
                format_func=format_strategy
            )
        else:
            selected_nuc_id = None
        
        # Bulk actions
        if nuc_ids:
            with st.expander("Bulk Actions", expanded=False):
                st.caption("Download multiple structures & data based on your current search.")
                
                download_mode = st.radio(
                    "Download Scope:",
                    ["Select Specific", "All Search Results"],
                    horizontal=True,
                    label_visibility="collapsed"
                )
                
                bulk_selected = []
                if download_mode == "Select Specific":
                    default_sel = nuc_ids[:5] if len(nuc_ids) > 0 else []
                    bulk_selected = st.multiselect(
                        "Select structures:",
                        nuc_ids,
                        default=default_sel,
                        format_func=format_strategy
                    )
                else:
                    bulk_selected = nuc_ids
               

