# ============================
# NucLigs Database – Streamlit App
# ============================

import streamlit as st
import py3Dmol
import pandas as pd
import io, zipfile, re

from supabase import create_client
from rdkit import Chem
from rdkit.Chem import (
    Descriptors, Crippen, rdMolDescriptors,
    Lipinski, QED
)

# ----------------------------------------------------
# PAGE CONFIG
# ----------------------------------------------------
st.set_page_config(
    page_title="NucLigs Database",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----------------------------------------------------
# CSS
# ----------------------------------------------------
st.markdown("""
<style>
.block-container { padding-top: 3.5rem; padding-bottom: 3rem; }

.meta-scroll { max-height: 70vh; overflow-y: auto; padding-right: 10px; }

.feature-card {
    background: #fff;
    border: 1px solid #e0e0e0;
    border-radius: 10px;
    padding: 15px;
    margin-bottom: 15px;
}

.ref-card {
    background: #fcfcfc;
    border-left: 4px solid #3498db;
    padding: 15px;
    margin-bottom: 15px;
    border-radius: 6px;
}

.data-row {
    display: flex;
    align-items: center;
    margin-bottom: 6px;
}

.data-label {
    min-width: 90px;
    font-size: 0.85rem;
suggest:
    font-weight: 600;
    color: #666;
}

.data-value {
    font-family: monospace;
    font-size: 0.9rem;
    color: #222;
}

.id-card {
    background: #f8f9fa;
    border-left: 5px solid #4CAF50;
    padding: 15px;
    border-radius: 8px;
    margin-bottom: 15px;
}

.badge-pass {
    background: #d4edda;
    color: #155724;
    padding: 4px 10px;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 700;
}

.badge-fail {
    background: #f8d7da;
    color: #721c24;
    padding: 4px 10px;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 700;
}
</style>
""", unsafe_allow_html=True)

# ----------------------------------------------------
# SUPABASE
# ----------------------------------------------------
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

# ----------------------------------------------------
# HELPERS
# ----------------------------------------------------
def format_pubmed(val):
    if not val or str(val).lower() == "nan":
        return ""
    m = re.search(r"\d+", str(val))
    return m.group(0) if m else ""

def ext_link(label, value, prefix):
    if not value:
        return
    v = str(value).strip()
    st.markdown(
        f"""
        <div class="data-row">
            <span class="data-label">{label}</span>
            <span class="data-value">
                <a href="{prefix}{v}" target="_blank">{v}</a>
            </span>
        </div>
        """,
        unsafe_allow_html=True
    )

def render_row(label, value):
    st.markdown(
        f"""
        <div class="data-row">
            <span class="data-label">{label}</span>
            <span class="data-value">{value}</span>
        </div>
        """,
        unsafe_allow_html=True
    )

# ----------------------------------------------------
# LOAD DATA
# ----------------------------------------------------
@st.cache_data
def load_metadata():
    b = supabase.storage.from_(BUCKET_META).download(META_FILE)
    df = pd.read_excel(io.BytesIO(b))
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    return df

@st.cache_data
def load_refs():
    b = supabase.storage.from_(BUCKET_META).download(REF_FILE)
    df = pd.read_excel(io.BytesIO(b))
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    return df

@st.cache_data
def load_pdb(name):
    b = supabase.storage.from_(BUCKET_PDB).download(name)
    return b.decode()

# ----------------------------------------------------
# 3D VIEWER WITH PNG BUTTON
# ----------------------------------------------------
def show_3d(pdb_text):
    v = py3Dmol.view(width=900, height=600)
    v.addModel(pdb_text, "pdb")
    v.setStyle({"stick": {"colorscheme": "greenCarbon"}})
    v.zoomTo()

    html = v._make_html().replace(
        "</body>",
        """
        <div style="text-align:center;margin-top:10px;">
          <button onclick="viewer.png()"
           style="padding:6px 14px;border-radius:6px;border:1px solid #ccc;">
           📷 Save PNG Snapshot
          </button>
        </div>
        </body>
        """
    )
    st.components.v1.html(html, height=650)

# ----------------------------------------------------
# MAIN APP
# ----------------------------------------------------
meta = load_metadata()
refs = load_refs()

st.sidebar.title("NucLigs Database")

nid = st.sidebar.selectbox(
    "Select Structure",
    meta["nl"].astype(str).tolist()
)

row = meta[meta["nl"].astype(str) == nid].iloc[0]
pdb_name = row["pdbs"]
pdb_text = load_pdb(pdb_name)

# ----------------------------------------------------
# LAYOUT
# ----------------------------------------------------
left, right = st.columns([1.4, 1])

with left:
    st.subheader(f"3D Structure – {nid}")
    show_3d(pdb_text)

with right:
    tabs = st.tabs(["Chemical Analysis", "Metadata", "References"])

    # ---------------- CHEMICAL ----------------
    with tabs[0]:
        mol = Chem.MolFromSmiles(row.get("smiles", "")) if pd.notna(row.get("smiles")) else None
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            viol = sum([
                mw > 500,
                logp > 5,
                Lipinski.NumHDonors(mol) > 5,
                Lipinski.NumHAcceptors(mol) > 10
            ])
            badge = "badge-pass" if viol == 0 else "badge-fail"

            st.markdown(f"<span class='{badge}'>Lipinski Violations: {viol}</span>", unsafe_allow_html=True)
            render_row("MW", f"{mw:.2f}")
            render_row("LogP", f"{logp:.2f}")
            render_row("TPSA", f"{rdMolDescriptors.CalcTPSA(mol):.2f}")
            render_row("QED", f"{QED.qed(mol):.3f}")

    # ---------------- METADATA ----------------
    with tabs[1]:
        st.markdown("<div class='id-card'>", unsafe_allow_html=True)
        render_row("NucL ID", nid)
        st.markdown("</div>", unsafe_allow_html=True)

        ext_link("ChEMBL", row.get("chembl_id"), "https://www.ebi.ac.uk/chembl/compound_report_card/")
        ext_link("DrugBank", row.get("drugbank_id"), "https://go.drugbank.com/drugs/")
        ext_link("PubMed", format_pubmed(row.get("pubmed_id")), "https://pubmed.ncbi.nlm.nih.gov/")

    # ---------------- REFERENCES ----------------
    with tabs[2]:
        chembl = str(row.get("chembl_id", "")).upper()
        pdbid = str(row.get("pdbs", "")).upper()

        hits = pd.DataFrame()
        if chembl and "chembl_id" in refs.columns:
            hits = refs[refs["chembl_id"].str.upper() == chembl]
        if hits.empty and "pdbs" in refs.columns:
            hits = refs[refs["pdbs"].str.upper() == pdbid]

        if hits.empty:
            st.info("No references found.")
        else:
            for _, r in hits.iterrows():
                st.markdown("<div class='ref-card'>", unsafe_allow_html=True)
                st.markdown(f"**{r.get('title','')}**")
                st.markdown(f"{r.get('journal','')} ({r.get('year','')})")
                pmid = format_pubmed(r.get("pubmed_id"))
                if pmid:
                    st.markdown(f"PubMed: [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid})")
                if pd.notna(r.get("doi")):
                    st.markdown(f"DOI: https://doi.org/{r['doi']}")
                st.markdown("</div>", unsafe_allow_html=True)
