import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(layout="wide", page_title="MOHEET v1.0 | PTP1B Hunter")

@st.cache_data
def load_data():
    try:
        # POINTING TO THE FAMILY TREE NOW
        df = pd.read_csv("data/processed/ptp1b_allosteric_family.csv")
        return df
    except FileNotFoundError:
        return None

df = load_data()

st.title("MOHEET v1.0: Cryptic Pocket Hunter")
st.markdown("**Target:** PTP1B (Allosteric) | **Series:** Benzofuran Family")

if df is not None:
    with st.sidebar:
        st.header("🎯 Filters")
        st.info(f"Family Members: {len(df)}")
        # Increased ranges so you don't filter out your own hits!
        mw_limit = st.slider("Max MW", 200, 800, 700) 
        logp_limit = st.slider("Max LogP", 0.0, 8.0, 7.0)
        
        mask = (df['mw'] <= mw_limit) & (df['logp'] <= logp_limit)
        filtered_df = df[mask]

    tab1, tab2 = st.tabs(["📊 Family Data", "🧪 Structure Gallery"])

    with tab1:
        st.dataframe(filtered_df, use_container_width=True, height=600)

    with tab2:
        st.write("### The Benzofuran Series")
        if not filtered_df.empty:
            cols = st.columns(3)
            for i, row in filtered_df.iterrows():
                with cols[i % 3]:
                    mol = Chem.MolFromSmiles(row['smiles'])
                    if mol:
                        img = Draw.MolToImage(mol, size=(250, 200))
                        st.image(img, caption=f"{row['chembl_id']}\nIC50: {row['ic50_nm']} nM\nSim: {round(row['similarity'],1)}%")
else:
    st.error("❌ Data not found. Run agents/cloner.py first.")