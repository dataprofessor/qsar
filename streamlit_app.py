import streamlit as st
import pandas as pd
from padelpy import padeldescriptor

st.title('💊 Bioactivity prediction app')

molecule = pd.read_csv('data/molecule.smi', header=None)
ic50 = pd.read_csv('data/hcv_ns5b_ic50_nm.csv', header=None)

example_input = 'CC(=O)OC1=CC=CC=C1C(=O)O'
smiles_txt = st.text_input('Enter SMILES notation', example_input)
st.info(smiles_txt)


fingerprint = 'Substructure'

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
st.write(descriptors)
