import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from bioinfokit import analys, visuz

# Read the DataFrame with p-value and log2 fold change columns
ucla_df = pd.read_csv("ucla_supplement_DEG_table_categorizedbycelltypes.csv")

# Convert the specified columns to numeric, errors='coerce' will replace non-numeric values with NaN
ucla_df['Log.fold.change'] = pd.to_numeric(ucla_df['Log.fold.change'], errors='coerce')
ucla_df['P Value'] = pd.to_numeric(ucla_df['P Value'], errors='coerce')

# Remove rows with NaN values in either of the two columns
ucla_df['Log10 P'] = ucla_df['P Value'].apply(lambda x: -np.log10(x))
ucla_df = ucla_df.replace([np.inf, -np.inf], np.nan)
ucla_df = ucla_df.dropna(subset=['Log.fold.change', 'P Value','Log10 P'])



# input
log = st.slider('Choose your log threshold', 0.0, 5.0, 2.5)
st.write("You selected: ", log)

p = st.slider('Choose your p threshold', 0.0, 0.2, 0.05)
st.write("You selected: ", p)


# Create a scatter plot
st.title("Volcano Plot")
plt.scatter(ucla_df["Log.fold.change"], ucla_df["Log10 P"], s=5)
plt.xlabel("Log2 Fold Change")
#plt.xlim([-5, 5])
plt.ylabel("P-Value")
#plt.ylim([0, 0.02])

# lines
plt.axvline(-log,color="grey",linestyle="--")
plt.axvline(log,color="grey",linestyle="--")
plt.axhline(-np.log10(p),color="grey",linestyle="--")

# find down- or up- regulated genes
down = ucla_df[(ucla_df['Log.fold.change']<=-log)&(ucla_df['P Value']<=p)]
up = ucla_df[(ucla_df['Log.fold.change']>=log)&(ucla_df['P Value']<=p)]

# highlight genes on plot
if len(down) >0:
    plt.scatter(x=down['Log.fold.change'],y=down["Log10 P"],s=8,label="Down-regulated",color="blue")
    # highlight genes on plot
if len(up) > 0:
    plt.scatter(x=up['Log.fold.change'],y=up["Log10 P"],s=8,label="Up-regulated",color="red")

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