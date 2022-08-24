import streamlit as st
import pandas as pd
from padelpy import padeldescriptor

st.title('💊 Bioactivity prediction app')

df = pd.read_csv('data/hcv_ns5b_curated_data.csv')
st.write(df)

example_input = "CC(=O)OC1=CC=CC=C1C(=O)O"
st.write(example_input)

smiles_txt = st.text_input('Enter SMILES notation', '')
st.info(smiles_txt)
      
