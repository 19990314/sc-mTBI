import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from bioinfokit import analys, visuz

# Read the DataFrame with p-value and log2 fold change columns
our_df = pd.read_csv("./DEG_day1_allEC_fc1.csv")
# print(our_df)
# Convert the specified columns to numeric, errors='coerce' will replace non-numeric values with NaN
our_df['avg_log2FC'] = pd.to_numeric(our_df['avg_log2FC'], errors='coerce')
our_df['p_val'] = pd.to_numeric(our_df['p_val'], errors='coerce')

# Remove rows with NaN values in either of the two columns
#our_df['Log10 P'] = our_df['p_val'].apply(lambda x: -np.log10(x))
our_df['Log10 P'] = our_df['p_val']
our_df = our_df.replace([np.inf, -np.inf], np.nan)
our_df = our_df.dropna(subset=['avg_log2FC', 'p_val', 'Log10 P'])

# input
log = 2.5
p = 0.05

p = st.slider('Choose your p threshold', 0.0, max(our_df['p_val']), 0.05)
st.write("You selected: ", p)

log = st.slider('Choose your (absolute) fold-change threshold', 0.0, max(our_df['avg_log2FC']), 2.5)
st.write("You selected: ", log)

gene_name = st.checkbox('Display Gene Names')

# realtime color
#colr = np.where((((our_df['Log10 P'] < p) & (our_df['avg_log2FC'] < -log)) | ((our_df['Log10 P'] < p) & (our_df['avg_log2FC'] > log))),"red", "grey")

# Create a scatter plot
st.title("Volcano Plot")
#plt.scatter(our_df["avg_log2FC"], our_df["Log10 P"], s=5, color = colr)
plt.scatter(our_df["avg_log2FC"], our_df["Log10 P"], s=5)
plt.xlabel("Log2 Fold Change")
# plt.xlim([-5, 5])
plt.ylabel("P-Value")
# plt.ylim([0, 0.02])

# lines
plt.axvline(-log, color="grey", linestyle="--")
plt.axvline(log, color="grey", linestyle="--")
plt.axhline(p, color="grey", linestyle="--")

# find down- or up- regulated genes
down = our_df[(our_df['avg_log2FC'] <= -log) & (our_df['Log10 P'] <= p)]
up = our_df[(our_df['avg_log2FC'] >= log) & (our_df['Log10 P'] <= p)]

# highlight genes on plot
if len(down) > 0:
    plt.scatter(x=down['avg_log2FC'], y=down["Log10 P"], s=8, label="Down-regulated", color="blue")
    # highlight genes on plot
if len(up) > 0:
    plt.scatter(x=up['avg_log2FC'], y=up["Log10 P"], s=8, label="Up-regulated", color="red")

# label names
if gene_name:
    for i, r in up.iterrows():
        plt.text(s=r['Gene'], x=r['avg_log2FC'], y=r["p_val"])
    for i, r in down.iterrows():
        plt.text(s=r['Gene'], x=r['avg_log2FC'], y=r["p_val"])

# final plot
st.pyplot(plt)

st.write("Up-regulated: ")
st.dataframe(up[['Gene', 'avg_log2FC', 'p_val']])

st.write("Down-regulated: ")
st.dataframe(down[['Gene', 'avg_log2FC', 'p_val']])
