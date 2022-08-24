import streamlit as st
import os
import pandas as pd
from padelpy import padeldescriptor

st.set_page_config(
  page_title='Bioactivity prediction app',
  page_icon='💊',
  layout='wide',
  initial_sidebar_state='collapsed')

if 'example_input' not in st.session_state:
  st.session_state.example_input = ''

if os.path.isfile('molecule.smi'):
  os.remove('molecule.smi') 
  
st.title('💊 Bioactivity prediction app')


# Input SMILES
st.sidebar.subheader('Input SMILES')

def insert_example_smiles():
    st.session_state.example_input = 'CC(=O)OC1=CC=CC=C1C(=O)O'
def clear_smiles():
    st.session_state.example_input = ''

smiles_txt = st.sidebar.text_input('Enter SMILES notation', st.session_state.example_input)
st.sidebar.button('Example input', on_click=insert_example_smiles)
st.sidebar.button('Clear input', on_click=clear_smiles)

if st.session_state.example_input == '':
  st.info('👈 Enter SMILES notation in the sidebar to proceed', icon='ℹ️')
else:
  st.info(smiles_txt, icon='ℹ️')

f = open('molecule.smi', 'w')
f.write(f'{smiles_txt}\tmol_001')
f.close()


# Compute PADEL descriptors
if os.path.isfile('molecule.smi'):
  padeldescriptor(mol_dir='molecule.smi', 
                  d_file='descriptors.csv',
                  descriptortypes='data/PubchemFingerprinter.xml', 
                  detectaromaticity=True,
                  standardizenitro=True,
                  standardizetautomers=True,
                  threads=2,
                  removesalt=True,
                  log=True,
                  fingerprints=True)

# Display PADEL descriptors
st.subheader('🔢 Descriptors')
descriptors = pd.read_csv('descriptors.csv')
descriptors.drop('Name', axis=1, inplace=True)
st.write(descriptors)
