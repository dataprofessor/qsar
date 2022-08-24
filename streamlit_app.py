import streamlit as st
import pandas as pd
from padelpy import padeldescriptor

if 'example_input' not in st.session_state:
  st.session_state.example_input = ''

st.title('💊 Bioactivity prediction app')

# Input SMILES
def insert_example_smiles():
    st.session_state.example_input = 'CC(=O)OC1=CC=CC=C1C(=O)O'
def clear_smiles():
    st.session_state.example_input = ''
    
smiles_txt = st.text_input('Enter SMILES notation', st.session_state.example_input)

button_col1, button_col2, button_col3 = st.columns(3)
button_col1.button('Insert example input', on_click=insert_example_smiles)
button_col2.button('Clear input', on_click=clear_url)

st.info(smiles_txt)

f = open('molecule.smi', 'w')
f.write(f'{Esmiles_txt}\t')


# Compute PADEL descriptors
padeldescriptor(mol_dir='data/molecule.smi', 
                d_file='data/descriptors.csv',
                descriptortypes='data/PubchemFingerprinter.xml', 
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=2,
                removesalt=True,
                log=True,
                fingerprints=True)

descriptors = pd.read_csv('descriptors.csv')
st.write(descriptors)
