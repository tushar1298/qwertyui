import streamlit as st
from supabase import create_client

# ==============================
# CONFIGURATION
# ==============================
SUPABASE_URL = "https://heuzgnhlrumyfcfigoon.supabase.co"

# ⚠️ MOVE THIS TO STREAMLIT SECRETS (see below)
SUPABASE_KEY = st.secrets["SUPABASE_KEY"]

BUCKET_NAME = "codes"

# ⚠️ UPDATE PATH IF FILE IS IN A FOLDER
SOURCE_FILENAME = "app1.py"   # e.g. "backend/app1.py"

# ==============================
# LOADER
# ==============================
def main():
    try:
        st.info("🔄 Loading application from Supabase...")

        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

        data = supabase.storage.from_(BUCKET_NAME).download(SOURCE_FILENAME)

        if not data:
            raise RuntimeError("Downloaded file is empty")

        source_code = data.decode("utf-8")

        # Execute remote app
        exec(source_code, globals())

    except Exception as e:
        st.error("❌ Failed to load the application from Supabase.")
        st.exception(e)

if __name__ == "__main__":
    main()
