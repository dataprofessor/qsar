import streamlit as st
import pandas as pd
from padelpy import padeldescriptor

st.title('💊 Bioactivity prediction app')

molecule = pd.read_csv('data/molecule.smi', header=None)
ic50 = pd.read_csv('data/hcv_ns5b_ic50_nm.csv', header=None)
st.write(molecule)
st.write(ic50)

example_input = "CC(=O)OC1=CC=CC=C1C(=O)O"
#st.write(example_input)

smiles_txt = st.text_input('Enter SMILES notation', '')
st.info(smiles_txt)
      
