import os
from collections import Counter
from typing import List

import pandas as pd
import streamlit as st
import plotly.express as px

import streamlit_permalink as stp

def get_env_str(var_name, default):
    return os.getenv(var_name, default)

PDB_APP_URL = get_env_str('PDB_APP_URL', 'https://pdb-cov.streamlit.app/')


import urllib.parse

def serialize_peptides(peptides: List[str]) -> str:
    """Serialize a list of peptides into a single string."""
    # counter
    peptide_counts = Counter(peptides)
    serialized = ','.join([f"{urllib.parse.quote(peptide, safe='')};{count}" for peptide, count in peptide_counts.items()])
    return serialized


with st.sidebar:
    parquet_file = st.file_uploader("Upload SagePro parquet file", type=['parquet', 'tsv'])

    if parquet_file:
        # Check file extension to determine how to read it
        file_extension = parquet_file.name.split('.')[-1].lower()
        
        if file_extension == 'parquet':
            df = pd.read_parquet(parquet_file)
        elif file_extension == 'tsv':
            df = pd.read_csv(parquet_file, sep='\t')
        else:
            st.error(f"Unsupported file type: {file_extension}")
            st.stop()
    else:
        st.stop()

    remove_decoys = st.checkbox("Remove decoys", value=True)

    qvalue_filter = st.multiselect("Q-value filter", options=['spectrum', 'peptide', 'protein'], default=['peptide'])
    qvalue_filter_value = st.number_input("Q-value filter value", value=0.01)

    keep_best_peptide = st.checkbox("Keep best peptide", value=False,
                                    help='Keep the best peptide, charge, filename pair')

    filter_protein = st.text_input("Filter Protein", value='', help='Filter proteins by substring match')

if remove_decoys:
    if 'is_decoy' in df:
        df = df[~df['is_decoy']]
    else:
        df = df[df['label']==1]

if len(qvalue_filter) > 0:
    if 'spectrum' in qvalue_filter:
        df = df[df['spectrum_q'] <= qvalue_filter_value]
    if 'peptide' in qvalue_filter:
        df = df[df['peptide_q'] <= qvalue_filter_value]
    if 'protein' in qvalue_filter:
        df = df[df['protein_q'] <= qvalue_filter_value]

if keep_best_peptide:
    # sort by hyperscore and keep the best peptide, charge, filename
    df = df.sort_values('hyperscore', ascending=False).groupby(['peptide', 'charge', 'filename']).head(1)

# calculate ppm from calcmass, expmass
df['ppm'] = (df['expmass'] - (df['calcmass'] + df['isotope_error'])) / df['expmass'] * 1e6

if filter_protein:
    # filter proteins by substring match
    df = df[df['proteins'].str.contains(filter_protein)]


c1, c2, c3 = st.columns([1, 1, 1])

# add metrics for number of peptides, proteins, and spectra
c1.metric(label="Peptides", value=len(df['peptide'].unique()))
c2.metric(label="Protein Groups", value=len(df['proteins'].unique()))
c3.metric(label="Spectra", value=len(df))


# plot a histogram of the peptide lengths (plotly)
tabs = st.tabs(['Data', 'Peptide Stats', 'Mass Error', 'Drift', 'Scatter', 'Coverage'])

deocy_col = 'is_decoy'
if 'is_decoy' not in df:
    deocy_col = 'label'

df['da'] = df['expmass'] - df['calcmass']

with tabs[0]:
    st.write(df)

with tabs[1]:
    fig = px.histogram(df, x='peptide_len', nbins=20, title='Peptide Lengths')
    st.plotly_chart(fig)

    #histogram of missed cleavages
    fig = px.histogram(df, x='missed_cleavages', title='Missed Cleavages')
    st.plotly_chart(fig)

    # charge
    fig = px.histogram(df, x='charge', title='Charge')
    st.plotly_chart(fig)

with tabs[2]:
    score_cos = ['spectrum_q', 'peptide_q', 'protein_q',
                 'sage_discriminant_score', 'poisson', 'delta_rt_model',
                 'hyperscore', 'matched_intensity_percent', 'delta_next']

    score_to_use = st.selectbox("Score to use", options=score_cos,
                                index=score_cos.index('hyperscore'))

    mass_error = st.selectbox('Mass Error Type', options=['ppm', 'da'], index=0)

    hover_data = ['peptide', 'charge', 'proteins']
    if 'ambiguity_sequence' in df.columns:
        hover_data.append('ambiguity_sequence')

    color = st.selectbox("Color by", options=[deocy_col, 'charge', 'missed_cleavages', 'peptide_len', 'proteins'], index=0)

    fig = px.scatter(df,
                     x=mass_error,
                     y=score_to_use,
                     title='Precursor Mass Error (ppm)',
                     color=color,
                     hover_data=hover_data)
    st.plotly_chart(fig)

with tabs[3]:
    # loop over files and plot drift
    for file in df['filename'].unique():
        file_df = df[df['filename'] == file]
        fig = px.scatter(file_df, x='rt', y='delta_rt_model', title=f'Drift vs RT Error {file}', color=deocy_col,
                         hover_data=['peptide'])
        st.plotly_chart(fig)

        # plot ppm vs rt
        fig = px.scatter(file_df, x='rt', y='ppm', title=f'PPM vs RT {file}', color=deocy_col,
                         hover_data=['peptide'])
        st.plotly_chart(fig)

with tabs[4]:
    #2d or 3d scatter plot
    scatter_type = st.radio("Scatter type", ['2D', '3D'], horizontal=True)

    axis_cols = ['spectrum_q', 'peptide_q', 'protein_q',
                 'sage_discriminant_score', 'poisson', 'delta_rt_model',
                 'hyperscore', 'matched_intensity_percent', 'delta_next', 'ppm']

    color_cols = [deocy_col, 'charge', 'missed_cleavages', 'peptide_len', 'semi_enzymatic', 'filename']

    color_axis = st.selectbox("Color axis", options=color_cols, index=color_cols.index(deocy_col))

    x_axis = st.selectbox("X-axis", options=axis_cols, index=axis_cols.index('hyperscore'))
    y_axis = st.selectbox("Y-axis", options=axis_cols, index=axis_cols.index('delta_next'))

    if scatter_type == '2D':
        fig = px.scatter(df, x=x_axis, y=y_axis, title=f'{x_axis} vs {y_axis}', color=color_axis,
                         hover_data=['peptide'])
        st.plotly_chart(fig)

    else:

        z_axis = st.selectbox("Z-axis", options=axis_cols, index=axis_cols.index('ppm'))
        fig = px.scatter_3d(df, x=x_axis, y=y_axis, z=z_axis, title=f'{x_axis} vs {y_axis} vs {z_axis}',
                            color=color_axis,
                            hover_data=['peptide'])
        # make bigger
        fig.update_layout(width=1000, height=1000)

        st.plotly_chart(fig)

with tabs[5]:
    # coverage plot

    rev_string = st.text_input("Reverse string", value='rev_')

    df['proforma'] = df.apply(lambda x: f"{x['peptide']}/{x['charge']}", axis=1)

    # group by protein and make alist of all peptides
    protein_df = df.groupby(['proteins']).agg({'proforma': list}).reset_index()
    protein_df['proteins'] = protein_df['proteins'].str.split(';')

    # explode the df so that each protein is on its own row (sep by ;)
    protein_df = protein_df.explode('proteins')

    # merge rows which have the same protein (combin peptide lists)
    protein_df = protein_df.groupby('proteins').agg({'proforma': 'sum'}).reset_index()

    protein_df['Locus Comps'] = protein_df['proteins'].str.split('|')
    # drop cols where there are not 3 values
    protein_df.reset_index(drop=True, inplace=True)

    protein_df['Database'] = protein_df['Locus Comps'].apply(lambda x: x[0] if len(x) == 3 else None)
    protein_df['Protein'] = protein_df['Locus Comps'].apply(lambda x: x[1] if len(x) == 3 else None)
    protein_df['Gene'] = protein_df['Locus Comps'].apply(lambda x: x[2] if len(x) == 3 else None)
    protein_df['Reverse'] = protein_df['proteins'].str.contains(rev_string)

    def make_link(protein_id, serialized_peptides, reverse, proteins, sequence_count, spectrum_count):

        
        
        if protein_id is None:
            input_type = 'Protein+Sequence'
            
            return f'{PDB_APP_URL}?input_type={input_type}&peptides={serialized_peptides}&reverse_protein={reverse}&protein_sequence={stp.EMPTY_STRING_URL_VALUE}&title={proteins}&subtitle=Seqs: {sequence_count} | Specs: {spectrum_count}&auto_spin=False&strip_mods=False&filter_unique=False'

        else:
            input_type = 'Protein+ID'
            return f'{PDB_APP_URL}?input_type={input_type}&protein_id={protein_id}&peptides={serialized_peptides}&reverse_protein={reverse}&auto_spin=False&strip_mods=False&filter_unique=False'

    protein_df['Sequence Count'] = protein_df['proforma'].apply(lambda x: len(set(x)))
    protein_df['Spectrum Count'] = protein_df['proforma'].apply(len)

    protein_df['SerializedPeptides'] = protein_df['proforma'].apply(serialize_peptides)
    protein_df['Link'] = protein_df.apply(lambda x: make_link(x['Protein'], x['SerializedPeptides'], x['Reverse'], x['proteins'], x['Sequence Count'],x['Spectrum Count']), axis=1)

    # drop locus comps, serialized peptides
    protein_df.drop(columns=['Locus Comps', 'SerializedPeptides'], inplace=True)

    protein_df = protein_df[['Link', 'proteins', 'Sequence Count', 'Spectrum Count']]

    # sort by sequence count
    protein_df = protein_df.sort_values('Sequence Count', ascending=False)

    st.dataframe(data=protein_df,
                 hide_index=True,
                 column_config={
                     'Sequence Count': st.column_config.NumberColumn(width="small"),
                     'Spectrum Count': st.column_config.NumberColumn(width="small"),
                     'Link': st.column_config.LinkColumn(display_text="PDB Viewer")
                 },
                 use_container_width=True)


