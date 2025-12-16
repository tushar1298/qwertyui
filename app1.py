import streamlit as st
import py3Dmol
import pandas as pd
import io
import zipfile
import re

from supabase import create_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, QED

# ====================================================
# PAGE CONFIG
# ====================================================
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ====================================================
# CSS
# ====================================================
st.markdown("""
<style>
.block-container { padding-top: 3rem; padding-bottom: 3rem; }

.meta-scroll { max-height: 70vh; overflow-y: auto; padding-right: 10px; }

.ref-card { background:#fcfcfc; border-left:4px solid #3498db;
            padding:14px; margin-bottom:14px; border-radius:6px; }

.ref-title { font-weight:700; font-size:1rem; margin-bottom:6px; }
.ref-meta { font-size:0.85rem; color:#555; margin-bottom:4px; }

.data-row { display:flex; align-items:center; margin-bottom:4px; }
.data-label { min-width:60px; font-weight:600; color:#666; }
.data-value { font-family:monospace; }

.id-card { background:#f8f9fa; border-left:5px solid #4CAF50;
           padding:14px; border-radius:8px; margin-bottom:14px; }

.stTabs [aria-selected="true"] { background:#e8f5e9 !important; }
</style>
""", unsafe_allow_html=True)

# ====================================================
# SUPABASE CONFIG
# ====================================================
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"
SUPABASE_KEY = "sb_secret_UuFsAopmAmHrdvHf6-mGBg_X0QNgMF5"

BUCKET_PDB = "NucLigs_PDBs"
BUCKET_META = "codes"

META_FILE = "NucLigs_metadata.xlsx"
REF_FILE = "references.xlsx"

@st.cache_resource
def init_supabase():
    return create_client(SUPABASE_URL, SUPABASE_KEY)

supabase = init_supabase()

# ====================================================
# LOADERS
# ====================================================
@st.cache_data
def load_metadata():
    try:
        b = supabase.storage.from_(BUCKET_META).download(META_FILE)
        df = pd.read_excel(io.BytesIO(b))
        df.columns = df.columns.str.lower().str.replace(" ", "_")
        return df
    except:
        return pd.DataFrame()

@st.cache_data
def load_references():
    try:
        b = supabase.storage.from_(BUCKET_META).download(REF_FILE)
        df = pd.read_excel(io.BytesIO(b))
        df.columns = df.columns.str.lower().str.replace(" ", "_")
        return df
    except:
        return pd.DataFrame()

@st.cache_data
def load_pdb(pdb_name: str):
    """SAFE Supabase PDB loader"""
    if not pdb_name or str(pdb_name).lower() == "nan":
        return None
    name = str(pdb_name).strip()
    if not name.lower().endswith(".pdb"):
        name += ".pdb"
    try:
        data = supabase.storage.from_(BUCKET_PDB).download(name)
        return data.decode("utf-8")
    except:
        return None

# ====================================================
# HELPERS
# ====================================================
def norm_chembl(x):
    if not x or str(x).lower() == "nan":
        return ""
    x = str(x).upper().replace("CHEMBL_", "CHEMBL")
    if not x.startswith("CHEMBL"):
        x = "CHEMBL" + x
    return x

def clean_pubmed(x):
    if not x or str(x).lower() == "nan":
        return ""
    m = re.search(r"\d+", str(x))
    return m.group(0) if m else ""

def link(label, url):
    return f"<a href='{url}' target='_blank'>{label}</a>"

# ====================================================
# 3D VIEWER WITH PNG EXPORT
# ====================================================
def show_3d(pdb_text, style="Stick"):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_text, "pdb")

    if style == "Stick":
        view.setStyle({"stick": {}})
    elif style == "Sphere":
        view.setStyle({"sphere": {}})
    else:
        view.setStyle({"line": {}})

    view.zoomTo()
    html = view._make_html()

    html = html.replace(
        "</body>",
        """
        <div style="text-align:center;margin-top:8px">
        <button onclick="viewer.png()"
        style="padding:6px 12px;border-radius:6px;border:1px solid #ccc;">
        📷 Save PNG Snapshot
        </button></div></body>
        """
    )

    st.components.v1.html(html, height=720)

# ====================================================
# MAIN APP
# ====================================================
df = load_metadata()
refs = load_references()

st.sidebar.title("NucLigs Browser")

if df.empty:
    st.error("Metadata not available")
    st.stop()

ids = sorted(df["nl"].astype(str).unique())
nid = st.sidebar.selectbox("Select NucL ID", ids)

row = df[df["nl"].astype(str) == nid].iloc[0]
pdb_text = load_pdb(row["pdbs"])

if not pdb_text:
    st.error("PDB file missing in Supabase")
    st.stop()

# ====================================================
# MAIN LAYOUT
# ====================================================
left, right = st.columns([1.5, 1])

with left:
    st.subheader(f"3D Structure: {nid}")
    style = st.selectbox("Style", ["Stick", "Sphere", "Line"])
    show_3d(pdb_text, style)

with right:
    tab1, tab2, tab3 = st.tabs(["Metadata", "Chemistry", "References"])

    # ---------------- METADATA ----------------
    with tab1:
        st.markdown(f"""
        <div class='id-card'>
        <b>Name:</b> {row.get("names","")}<br><br>
        <b>DrugBank:</b> {link(row.get("drugbank_id","N/A"),
        f"https://go.drugbank.com/drugs/{row.get('drugbank_id')}")}<br>
        <b>ChEMBL:</b> {link(norm_chembl(row.get("chembl_id")),
        f"https://www.ebi.ac.uk/chembl/compound_report_card/{norm_chembl(row.get('chembl_id'))}")}<br>
        <b>PubMed:</b> {link(clean_pubmed(row.get("pubmed_id")),
        f"https://pubmed.ncbi.nlm.nih.gov/{clean_pubmed(row.get('pubmed_id'))}/")}
        </div>
        """, unsafe_allow_html=True)

    # ---------------- CHEMISTRY ----------------
    with tab2:
        mol = Chem.MolFromSmiles(str(row.get("smiles","")))
        if mol:
            st.write("Molecular Weight:", Descriptors.MolWt(mol))
            st.write("LogP:", Crippen.MolLogP(mol))
            st.write("TPSA:", rdMolDescriptors.CalcTPSA(mol))

    # ---------------- REFERENCES ----------------
    with tab3:
        chembl = norm_chembl(row.get("chembl_id"))
        pdb = str(row.get("pdbs")).replace(".pdb","")

        hits = pd.DataFrame()
        if not refs.empty:
            hits = refs[
                (refs.get("chembl_id","").apply(norm_chembl) == chembl) |
                (refs.get("pdbs","").astype(str) == pdb)
            ]

        if hits.empty:
            st.info("No references found")
        else:
            for _, r in hits.iterrows():
                pm = clean_pubmed(r.get("pubmed_id"))
                st.markdown(f"""
                <div class='ref-card'>
                <div class='ref-title'>{r.get("title","")}</div>
                <div class='ref-meta'>{r.get("journal","")} ({r.get("year","")})</div>
                <div class='data-row'><span class='data-label'>PubMed</span>
                <span class='data-value'>
                <a href='https://pubmed.ncbi.nlm.nih.gov/{pm}/' target='_blank'>{pm}</a>
                </span></div>
                <div class='data-row'><span class='data-label'>DOI</span>
                <span class='data-value'>
                <a href='https://doi.org/{r.get("doi","")}' target='_blank'>{r.get("doi","")}</a>
                </span></div>
                </div>
                """, unsafe_allow_html=True)
