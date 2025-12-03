import streamlit as st
import pandas as pd
import psycopg2
import requests
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol

# ---------------------------------------
# GitHub PDB SETTINGS
# ---------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/"

# ---------------------------------------
# Load PDB file list from GitHub Repo
# ---------------------------------------
@st.cache_data
def list_pdb_files():
    r = requests.get(GITHUB_API_URL)
    files = r.json()
    pdb_files = [f["name"] for f in files if f["name"].endswith(".pdb")]
    return pdb_files

# ---------------------------------------
# Fetch PDB File from GitHub
# ---------------------------------------
def fetch_pdb_from_github(filename: str) -> str | None:
    url = f"{GITHUB_RAW_BASE}{filename}"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            return r.text
        else:
            st.error("❌ Could not fetch PDB file.")
            return None
    except Exception as e:
        st.error(f"Error: {e}")
        return None

# ---------------------------------------
# Display 3D Structure
# ---------------------------------------
def show_3d_pdb(pdb_text: str):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=520)

# ---------------------------------------
# Database Connection (for molecule table)
# ---------------------------------------
@st.cache_resource
def get_connection():
    db = st.secrets["db"]
    conn = psycopg2.connect(
        host=db["host"],
        database=db["dbname"],
        user=db["user"],
        password=db["password"],
        port=db["port"],
    )
    return conn

@st.cache_data
def load_molecules():
    conn = get_connection()
    df = pd.read_sql("SELECT * FROM molecules;", conn)
    return df

# ---------------------------------------
# 2D Structure
# ---------------------------------------
def show_2d_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img, caption="2D Structure")
    else:
        st.warning("Invalid SMILES")

# ---------------------------------------
# STREAMLIT UI
# ---------------------------------------
st.title("🧬 Molecule Viewer (2D + 3D)")

# Load molecule data
df = load_molecules()

# -----------------------------
# Molecule selection
# -----------------------------
with st.sidebar:
    st.header("🔎 Search Molecule")
    text_query = st.text_input("Search by name or ID")

if text_query:
    mask = df["name"].astype(str).str.contains(text_query, case=False, na=False)
    results = df[mask]
else:
    results = df

if len(results) == 0:
    st.warning("No molecule found.")
else:
    selected = st.sidebar.selectbox(
        "Select molecule",
        results.index,
        format_func=lambda i: f"{results.at[i, 'molecule_id']} - {results.at[i, 'name']}"
    )

    row = results.loc[selected]

    st.subheader(f"🧪 {row['name']} ({row['molecule_id']})")
    st.markdown("### 📊 Properties")
    st.write(row)

    # -------------------------------
    # 2D Structure
    # -------------------------------
    st.markdown("### 🔷 2D Structure")
    show_2d_structure(row["smiles"])

    # -------------------------------
    # 3D PDB Viewer (from GitHub)
    # -------------------------------
    st.markdown("### 🧬 3D Structure (PDB from GitHub)")

    pdb_files = list_pdb_files()
    pdb_choice = st.selectbox("Select a PDB file", pdb_files)

    if pdb_choice:
        pdb_text = fetch_pdb_from_github(pdb_choice)
        if pdb_text:
            show_3d_pdb(pdb_text)
