import streamlit as st
import os
import pickle
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from stmol import speck_plot
from padelpy import padeldescriptor

st.set_page_config(
  page_title='BioAct',
  page_icon='💊',
  layout='wide',
  initial_sidebar_state='collapsed')

if 'smiles_input' not in st.session_state:
  st.session_state.smiles_input = ''

if os.path.isfile('molecule.smi'):
  os.remove('molecule.smi') 
  
st.title('💊 BioAct: A BioActivity prediction app')

with st.expander('About this app'):
  st.write('''
    This Biological Activity (**BioAct**) prediction app allow users to easily evaluate the putative biological activity of a query molecule against the target protein being investigated.
    
    This app is based on the following Python libraries:
    - `streamlit`
    - `pandas`
    - `rdkit`
    - `stmol`
    - `padelpy`
    ''')

# Input SMILES
st.sidebar.subheader('Input molecule')

def insert_example_smiles():
    st.session_state.smiles_input = 'CC(=O)OC1=CC=CC=C1C(=O)O'
def clear_smiles():
    st.session_state.smiles_input = ''

smiles_txt = st.sidebar.text_input('Enter SMILES notation', st.session_state.smiles_input)
st.sidebar.button('Example input', on_click=insert_example_smiles)
st.sidebar.button('Clear input', on_click=clear_smiles)


if st.session_state.smiles_input == '':
  st.subheader('Welcome to the app!')
  st.info('Enter SMILES notation in the sidebar to proceed', icon='👈')
else:
  st.subheader('⚛️ Input molecule:')
  st.info(smiles_txt)

f = open('molecule.smi', 'w')
f.write(f'{smiles_txt}\tmol_001')
f.close()

# Show molecule
if st.session_state.smiles_input != '':
  m = Chem.MolFromSmiles(smiles_txt)
  m2 = Chem.AddHs(m)
  AllChem.EmbedMolecule(m2,randomSeed=0xf00d)
  AllChem.MMFFOptimizeMolecule(m2)
  m3 = Chem.MolToXYZBlock(m2)
  Chem.MolToXYZFile(m2, 'molecule.xyz')

  mol_col1, mol_col2 = st.columns(2)
  
  with mol_col2:
    with st.expander('Show XYZ file content'):
      st.code(m3)

  with mol_col1:
    with st.expander('Show molecular structure via stmol'):
      f2 = open('molecule.xyz', 'r')
      molecule_xyz = f2.read()
      speck_plot(molecule_xyz, wbox_width='350px', wbox_height='600px')


# Compute PADEL descriptors
if st.session_state.smiles_input != '':
  st.subheader('🔢 Descriptors')
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

  descriptors = pd.read_csv('descriptors.csv')
  descriptors.drop('Name', axis=1, inplace=True)
  st.write(descriptors)
  
 
# Read in saved classification model
# if st.session_state.smiles_input != '':
  # model = pickle.load(open('data/oversampling_PubChem_RandomForestClassifier.pkl', 'rb'))
  # importances = model.feature_importances_
  # model_importances = pd.Series(importances, index=feature_names)
  # st.write(model_importances)
