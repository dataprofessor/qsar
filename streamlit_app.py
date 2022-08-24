import streamlit as st
import pandas as pd

st.title('💊 Bioactivity prediction app')

df = pd.read_csv('data/hcv_ns5b_curated_data.csv')

st.write(df)
