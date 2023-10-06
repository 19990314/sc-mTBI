import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from bioinfokit import analys, visuz

# Read the DataFrame with p-value and log2 fold change columns
our_df = pd.read_csv("../../mTBI_new/1006_DEG_allEC_0v1_p01_wilcox_ourdata.CSV")

# Convert the specified columns to numeric, errors='coerce' will replace non-numeric values with NaN
our_df['avg_log2FC'] = pd.to_numeric(our_df['avg_log2FC'], errors='coerce')
our_df['p_val'] = pd.to_numeric(our_df['p_val'], errors='coerce')

# Remove rows with NaN values in either of the two columns
our_df['Log10 P'] = our_df['p_val'].apply(lambda x: -np.log10(x))
our_df = our_df.replace([np.inf, -np.inf], np.nan)
our_df = our_df.dropna(subset=['avg_log2FC', 'p_val','Log10 P'])


# input
log = 0.02
p = 0.05

p = st.slider('Choose your p threshold', 0.0, 0.2, 0.05)
st.write("You selected: ", p)


# Create a scatter plot
st.title("Volcano Plot")
plt.scatter(our_df["avg_log2FC"], our_df["Log10 P"], s=5)
plt.xlabel("Log2 Fold Change")
#plt.xlim([-5, 5])
plt.ylabel("P-Value")
#plt.ylim([0, 0.02])

# lines
plt.axvline(-log,color="grey",linestyle="--")
plt.axvline(log,color="grey",linestyle="--")
plt.axhline(-np.log10(p),color="grey",linestyle="--")

# find down- or up- regulated genes
down = our_df[(our_df['avg_log2FC']<=-log)&(our_df['P Value']<=p)]
up = our_df[(our_df['avg_log2FC']>=log)&(our_df['P Value']<=p)]

# highlight genes on plot
if len(down) >0:
    plt.scatter(x=down['avg_log2FC'],y=down["Log10 P"],s=8,label="Down-regulated",color="blue")
    # highlight genes on plot
if len(up) > 0:
    plt.scatter(x=up['avg_log2FC'],y=up["Log10 P"],s=8,label="Up-regulated",color="red")

# label names
for i,r in up.iterrows():
    plt.text(s=r['Gene'], x=r['Log.fold.change'], y=r["Log10 P"])
for i,r in down.iterrows():
    plt.text(s=r['Gene'], x=r['Log.fold.change'], y=r["Log10 P"])

# final plot
st.pyplot(plt)


st.write("Up-regulated: ")
st.dataframe(up[['Gene', 'Log.fold.change', 'P Value']])

st.write("Down-regulated: ")
st.dataframe(down[['Gene', 'Log.fold.change', 'P Value']])