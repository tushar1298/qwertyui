import streamlit as st
import requests
import py3Dmol
import pandas as pd

# ----------------------------------------------------
# GitHub repo settings
# ----------------------------------------------------
GITHUB_API_URL = "https://api.github.com/repos/tushar1298/qwertyui/contents"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/"

# CHANGE THIS if your metadata filename is different
METADATA_URL = "https://raw.githubusercontent.com/tushar1298/qwertyui/main/NucLigs_data_2811.xlsx"
# If you use CSV instead, e.g. metadata.csv, just change the extension and it will work.


# ----------------------------------------------------
# Helper: list all PDB files in the GitHub repo
# ----------------------------------------------------
@st.cache_data
def list_pdb_files():
    r = requests.get(GITHUB_API_URL)
    r.raise_for_status()
    files = r.json()
    pdb_files = [f["name"] for f in files if isinstance(f, dict) and f.get("name", "").endswith(".pdb")]
    return sorted(pdb_files)


# ----------------------------------------------------
# Helper: fetch a PDB file from GitHub
# ----------------------------------------------------
def fetch_pdb_from_github(filename: str) -> str | None:
    url = f"{GITHUB_RAW_BASE}{filename}"
    try:
        r = requests.get(url)
        if r.status_code == 200 and r.text.strip():
            return r.text
        else:
            st.error(f"❌ Could not fetch PDB file: {filename}")
            return None
    except Exception as e:
        st.error(f"Error fetching PDB: {e}")
        return None


# ----------------------------------------------------
# Helper: show PDB in 3D stick format
# ----------------------------------------------------
def show_3d_pdb(pdb_text: str):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"stick": {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=520)


# ----------------------------------------------------
# Load metadata from Excel / CSV in GitHub
# ----------------------------------------------------
@st.cache_data
def load_metadata():
    # Decide based on file extension
    if METADATA_URL.endswith(".csv"):
        df = pd.read_csv(METADATA_URL)
    else:
        df = pd.read_excel(METADATA_URL)
    # Normalize column names to lowercase for easier matching
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def find_metadata_for_pdb(df: pd.DataFrame, pdb_filename: str) -> pd.DataFrame | None:
    """
    Try to match the selected pdb file with a row in metadata.

    We try these possibilities:
      - column 'pdb_file' equals 'name.pdb'
      - column 'file_name' / 'filename' equals 'name.pdb'
      - column 'pdb' or 'pdb_id' equals 'name' (without .pdb)
    """
    fname = pdb_filename.strip()
    root = fname[:-4] if fname.lower().endswith(".pdb") else fname

    candidates = []

    # full filename matches
    for col in ["pdb_file", "file_name", "filename"]:
        if col in df.columns:
            mask = df[col].astype(str).str.lower() == fname.lower()
            candidates.append(df[mask])

    # id (without .pdb) matches
    for col in ["pdb", "pdb_id"]:
        if col in df.columns:
            mask = df[col].astype(str).str.lower() == root.lower()
            candidates.append(df[mask])

    # pick the first non-empty candidate
    for cand in candidates:
        if cand is not None and not cand.empty:
            return cand

    return None


# ----------------------------------------------------
# STREAMLIT APP
# ----------------------------------------------------
st.title("🧬 GitHub PDB 3D Viewer with Metadata")

# Load list of PDBs and metadata
try:
    pdb_files = list_pdb_files()
except Exception as e:
    st.error(f"Error listing PDB files: {e}")
    pdb_files = []

try:
    metadata_df = load_metadata()
except Exception as e:
    st.error(f"Error loading metadata file: {e}")
    metadata_df = pd.DataFrame()

if not pdb_files:
    st.warning("No PDB files found in the GitHub repository.")
else:
    # Select PDB file
    selected_pdb = st.selectbox("Select a PDB file", pdb_files)

    if selected_pdb:
        pdb_text = fetch_pdb_from_github(selected_pdb)
        if pdb_text:
            # Two columns: left = 3D, right = metadata table
            col1, col2 = st.columns([2, 1])

            with col1:
                st.subheader(f"3D Structure: {selected_pdb}")
                show_3d_pdb(pdb_text)

            with col2:
                st.subheader("Metadata")
                if metadata_df.empty:
                    st.info("No metadata file loaded.")
                else:
                    matched = find_metadata_for_pdb(metadata_df, selected_pdb)
                    if matched is None or matched.empty:
                        st.info("No metadata found for this PDB.")
                    else:
                        # reset index and show as table
                        st.dataframe(matched.reset_index(drop=True))
