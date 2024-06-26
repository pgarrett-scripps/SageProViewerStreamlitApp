import os
from collections import Counter

import pandas as pd
import streamlit as st
import plotly.express as px

def get_env_str(var_name, default):
    return os.getenv(var_name, default)

PDB_APP_URL = get_env_str('PDB_APP_URL', 'http://pdb-coverage.streamlit.app/')

with st.sidebar:
    parquet_file = st.file_uploader("Upload SagePro parquet file", type=['parquet'])

    if not parquet_file:
        st.stop()

    df = pd.read_parquet(parquet_file)

    remove_decoys = st.checkbox("Remove decoys", value=True)

    qvalue_filter = st.multiselect("Q-value filter", options=['spectrum', 'peptide', 'protein'], default=['spectrum'])
    qvalue_filter_value = st.number_input("Q-value filter value", value=0.05)

    keep_best_peptide = st.checkbox("Keep best peptide", value=True,
                                    help='Keep the best peptide, charge, filename pair')

if remove_decoys:
    df = df[~df['is_decoy']]

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
df['ppm'] = (df['calcmass'] + df['isotope_error'] - df['expmass']) / df['calcmass'] * 1e6

# plot a histogram of the peptide lengths (plotly)
tabs = st.tabs(['Data', 'Peptide Stats', 'Mass Error', 'Drift', 'Scatter', 'Coverage'])

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

    fig = px.scatter(df,
                     x='ppm',
                     y=score_to_use,
                     title='Precursor Mass Error (ppm)',
                     color='is_decoy',
                     hover_data=['peptide'])
    st.plotly_chart(fig)

with tabs[3]:
    # loop over files and plot drift
    for file in df['filename'].unique():
        file_df = df[df['filename'] == file]
        fig = px.scatter(file_df, x='rt', y='delta_rt_model', title=f'Drift vs RT Error {file}', color='is_decoy',
                         hover_data=['peptide'])
        st.plotly_chart(fig)

        # plot ppm vs rt
        fig = px.scatter(file_df, x='rt', y='ppm', title=f'PPM vs RT {file}', color='is_decoy',
                         hover_data=['peptide'])
        st.plotly_chart(fig)

with tabs[4]:
    #2d or 3d scatter plot
    scatter_type = st.radio("Scatter type", ['2D', '3D'], horizontal=True)

    axis_cols = ['spectrum_q', 'peptide_q', 'protein_q',
                 'sage_discriminant_score', 'poisson', 'delta_rt_model',
                 'hyperscore', 'matched_intensity_percent', 'delta_next', 'ppm']

    color_cols = ['is_decoy', 'charge', 'missed_cleavages', 'peptide_len', 'semi_enzymatic', 'filename']

    color_axis = st.selectbox("Color axis", options=color_cols, index=color_cols.index('is_decoy'))

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

    # group by protein and make alist of all peptides
    protein_df = df.groupby(['proteins']).agg({'peptide': list}).reset_index()
    protein_df['proteins'] = protein_df['proteins'].str.split(';')

    # explode the df so that each protein is on its own row (sep by ;)
    protein_df = protein_df.explode('proteins')

    # merge rows which have the same protein (combin peptide lists)
    protein_df = protein_df.groupby('proteins').agg({'peptide': 'sum'}).reset_index()

    protein_df['Locus Comps'] = protein_df['proteins'].str.split('|')
    # drop cols where there are not 3 values
    protein_df = protein_df[protein_df['Locus Comps'].apply(len) == 3]
    protein_df.reset_index(drop=True, inplace=True)

    protein_df['Database'] = protein_df['Locus Comps'].apply(lambda x: x[0])
    protein_df['Protein'] = protein_df['Locus Comps'].apply(lambda x: x[1])
    protein_df['Gene'] = protein_df['Locus Comps'].apply(lambda x: x[2])
    protein_df['Reverse'] = protein_df['Database'].str.contains(rev_string)

    # serialize_peptides
    def serialize_peptides(peptides):
        peptide_counts = Counter(peptides)
        return ','.join([f'{k};{v}' for k, v in peptide_counts.items()])

    def make_link(protein_id, serialized_peptides, reverse):
        return f'{PDB_APP_URL}?protein_id={protein_id}&input={serialized_peptides}&input_type=redundant_peptides&reverse_protein={reverse}'

    protein_df['SerializedPeptides'] = protein_df['peptide'].apply(serialize_peptides)
    protein_df['Link'] = protein_df.apply(lambda x: make_link(x['Protein'], x['SerializedPeptides'], x['Reverse']), axis=1)

    # drop locus comps, serialized peptides
    protein_df.drop(columns=['Locus Comps', 'SerializedPeptides'], inplace=True)

    protein_df['Sequence Count'] = protein_df['peptide'].apply(lambda x: len(set(x)))
    protein_df['Spectrum Count'] = protein_df['peptide'].apply(len)

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


