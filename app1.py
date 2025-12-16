import streamlit as st
import py3Dmol
import pandas as pd
import io, zipfile, re
from supabase import create_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED

# ----------------------------------------------------
# Page setup
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# CSS Styling
# ----------------------------------------------------
st.markdown("""<style>
.block-container{padding-top:3.5rem;padding-bottom:3rem}
.meta-scroll{max-height:70vh;overflow-y:auto;padding-right:12px}
.feature-card{background:#fff;border:1px solid #e0e0e0;border-radius:10px;padding:15px;margin-bottom:15px}
.data-row{display:flex;justify-content:space-between;margin-bottom:8px}
.data-label{font-size:.85rem;color:#666;font-weight:600}
.data-value{font-family:monospace;font-size:.9rem;font-weight:600}
.id-card{background:#f8f9fa;border-left:5px solid #4CAF50;padding:15px;border-radius:8px;margin-bottom:15px}
</style>""", unsafe_allow_html=True)

# ----------------------------------------------------
# Supabase Configuration
# ----------------------------------------------------
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_publishable_AM951Hs4gISMnct_hoTOkA_CnjMPj97"

BUCKET_NAME = "NucLigs_PDBs"
METADATA_BUCKET = "codes"
METADATA_FILENAME = "NucLigs_metadata.xlsx"
METADATA_REF_FILENAME = "references.xlsx"

@st.cache_resource
def init_supabase():
    return create_client(SUPABASE_URL, SUPABASE_KEY)

supabase = init_supabase()

# ----------------------------------------------------
# Helpers
# ----------------------------------------------------
def render_row(label, value):
    st.markdown(
        f"<div class='data-row'><span class='data-label'>{label}</span>"
        f"<span class='data-value'>{value}</span></div>",
        unsafe_allow_html=True
    )

def render_link_row(label, value, base_url):
    if not value or str(value).lower() == "nan":
        return
    st.markdown(
        f"<div class='data-row'><span class='data-label'>{label}</span>"
        f"<span class='data-value'><a href='{base_url}{value}' target='_blank'>{value}</a></span></div>",
        unsafe_allow_html=True
    )

def fetch_pdb(fname):
    return supabase.storage.from_(BUCKET_NAME).download(fname).decode()

@st.cache_data
def load_metadata():
    data = supabase.storage.from_(METADATA_BUCKET).download(METADATA_FILENAME)
    df = pd.read_excel(io.BytesIO(data))
    df.columns = df.columns.str.lower()
    return df

@st.cache_data
def load_refs():
    data = supabase.storage.from_(METADATA_BUCKET).download(METADATA_REF_FILENAME)
    return pd.read_excel(io.BytesIO(data))

# ----------------------------------------------------
# Main App
# ----------------------------------------------------
df = load_metadata()
refs_df = load_refs()

selected = st.sidebar.selectbox("Select NucLigs ID", df["nl"].astype(str))
row = df[df["nl"].astype(str) == selected].iloc[0]
data = row.to_dict()
pdb_text = fetch_pdb(row["pdbs"])

col_left, col_right = st.columns([1.5, 1])

# LEFT
with col_left:
    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    st.components.v1.html(view._make_html(), height=680)

# RIGHT
with col_right:
    tab_analysis, tab_metadata, tab_refs = st.tabs(
        ["Chemical Analysis", "Metadata Record", "References"]
    )

    # ---------------- Metadata Tab ----------------
    with tab_metadata:
        st.markdown("<div class='meta-scroll'>", unsafe_allow_html=True)

        st.markdown(
            f"<div class='id-card'><b>NucLigs ID:</b> {data.get('nl')}<br>"
            f"{data.get('names','')}</div>",
            unsafe_allow_html=True
        )

        st.markdown("<div class='feature-card'><h5>General Info</h5>", unsafe_allow_html=True)

        render_link_row("PubChem ID", data.get("pubchem_id"),
                        "https://pubchem.ncbi.nlm.nih.gov/compound/")
        render_link_row("DrugBank ID", data.get("drugbank_id"),
                        "https://go.drugbank.com/drugs/")
        render_link_row("ChEMBL ID", data.get("chembl_id"),
                        "https://www.ebi.ac.uk/chembl/compound_report_card/")

        skip = {"nl","names","pdbs","pubchem_id","drugbank_id","chembl_id"}
        for k, v in data.items():
            if k not in skip and pd.notna(v):
                render_row(k.replace("_"," ").title(), v)

        st.markdown("</div></div>", unsafe_allow_html=True)

    # ---------------- References Tab ----------------
    with tab_refs:
        st.markdown("<div class='meta-scroll'>", unsafe_allow_html=True)
        pdb_id = str(data.get("pdbs"))

        matches = refs_df[refs_df["pdbs"].astype(str) == pdb_id] if not refs_df.empty else []

        if len(matches):
            for _, r in matches.iterrows():
                st.markdown("<div class='feature-card'>", unsafe_allow_html=True)
                render_row("Title", r.get("title",""))
                if pd.notna(r.get("doi")):
                    st.markdown(
                        f"<a href='https://doi.org/{r.get('doi')}' target='_blank'>DOI</a>",
                        unsafe_allow_html=True
                    )
                st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.info("No references found.")

        st.markdown("</div>", unsafe_allow_html=True)

st.success("✔ PubChem, DrugBank & ChEMBL IDs are clickable")
