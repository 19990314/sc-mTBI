import os
import pandas as pd

deg_tracker = {
    'studies': ["Sham v.s. Day1", "Sham v.s. Day3", "Sham v.s. Day7"],
    'upregulated': [],
    'downregulated': []
}

def get_EC_celltype(filename):
    for i in filename.split("_"):
        if "EC" in i:
            return i

def read_deg_tables(filename_prefix):
    # Directory containing the CSV files
    directory = '../../mTBI_scRNA_seq/DEGs/'

    # Get a list of CSV files with the specified prefix
    csv_files = [file for file in os.listdir(directory) if file.startswith(filename_prefix)]

    # Create an empty dictionary to store dataframes
    dataframes = {}

    # Loop through the CSV files and read them into dataframes
    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        # Use the file name (without extension) as the dataframe key
        dataframe_key = os.path.splitext(csv_file)[0]
        # Read the CSV file into a dataframe
        dataframes[dataframe_key] = pd.read_csv(file_path, header=0, names = ['gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj'], dtype='object')
        dataframes[dataframe_key].set_index('gene', inplace=True)
        dataframes[dataframe_key] = dataframes[dataframe_key].apply(pd.to_numeric, errors='coerce')

    return dataframes


def filter_degs(dataframes, pval, fc):
    for df_key in dataframes:
        size_tracker['celltype'].append(get_EC_celltype(df_key))
        # Filter rows/degs
        size_tracker['before filtering'].append(len(dataframes[df_key]))
        dataframes[df_key] = dataframes[df_key][dataframes[df_key]['p_val'] <= pval]
        size_tracker['after p-val'].append(len(dataframes[df_key]))
        dataframes[df_key] = dataframes[df_key][abs(dataframes[df_key]['avg_log2FC']) >= fc]
        size_tracker['after fc'].append(len(dataframes[df_key]))
        if "Tmem252" in dataframes[df_key].index:
            size_tracker['Tmem252 check'].append("Yes")
        else:
            size_tracker['Tmem252 check'].append("No")
    return dataframes


def find_common_degs(dataframes):
    common_degs = ()
    for df_key in dataframes:
        if len(common_degs) == 0:
            common_degs = set(dataframes[df_key].index)
            continue
        common_degs = set(common_degs.intersection(set(dataframes[df_key].index)))

    return list(common_degs)


def diff_up_vs_downregulation(dataframes, common_degs):
    up = []
    down = []
    for gene in common_degs:
        concordance_flag = 0
        regulation = ""
        for df_key in dataframes:
             if dataframes[df_key].loc[gene]['avg_log2FC'] < 0:
                 if regulation == "":
                     regulation = "down"
                 elif regulation == "down":
                     continue
                 elif regulation == "up":
                     concordance_flag = 1
                     break
             elif dataframes[df_key].loc[gene]['avg_log2FC'] > 0:
                 if regulation == "":
                     regulation = "up"
                 elif regulation == "up":
                     continue
                 elif regulation == "down":
                     concordance_flag = 1
                     break
        if concordance_flag == 0:
            if regulation == "up":
                up.append(gene)
            else:
                down.append(gene)
    return up, down



## ---------day1---------
# DEG size tracker
size_tracker = {
        'celltype': [],
        'before filtering': [],
        'after p-val': [],
        'after fc': [],
        'Tmem252 check': []
}

# read and filter
degs_for_day1 = read_deg_tables("1006_DEG_0v1")
degs_for_day1 = filter_degs(degs_for_day1, 0.001, 0.05)

# diff degs
common_degs = find_common_degs(degs_for_day1)
#print(common_degs)
up, down = diff_up_vs_downregulation(degs_for_day1, common_degs)
deg_tracker['upregulated'].append(up)
deg_tracker['downregulated'].append(down)

# output
df = pd.DataFrame(size_tracker)
df.to_csv('../../mTBI_scRNA_seq/DEGs/day1_deg_counts.csv', index=False)





## ---------day3---------
# DEG size tracker
size_tracker = {
        'celltype': [],
        'before filtering': [],
        'after p-val': [],
        'after fc': [],
        'Tmem252 check': []
}

# read and filter
degs_for_day3 = read_deg_tables("1005_DEG_0v3")
degs_for_day3 = filter_degs(degs_for_day3, 0.001, 0.05)

# diff degs
common_degs = find_common_degs(degs_for_day3)
#print(common_degs)
up, down = diff_up_vs_downregulation(degs_for_day3, common_degs)
deg_tracker['upregulated'].append(up)
deg_tracker['downregulated'].append(down)

# output
df = pd.DataFrame(size_tracker)
df.to_csv('../../mTBI_scRNA_seq/DEGs/day3_deg_counts.csv', index=False)




## ---------day7---------
# DEG size tracker
size_tracker = {
        'celltype': [],
        'before filtering': [],
        'after p-val': [],
        'after fc': [],
        'Tmem252 check': []
}

# read and filter
degs_for_day7 = read_deg_tables("1005_DEG_0v7")
degs_for_day7 = filter_degs(degs_for_day7, 0.001, 0.05)

# diff degs
common_degs = find_common_degs(degs_for_day7)
#print(common_degs)
up, down = diff_up_vs_downregulation(degs_for_day7, common_degs)
deg_tracker['upregulated'].append(up)
deg_tracker['downregulated'].append(down)

# output
df = pd.DataFrame(size_tracker)
df.to_csv('../../mTBI_scRNA_seq/DEGs/day7_deg_counts.csv', index=False)



# output degs for day 1,3,7
df = pd.DataFrame(deg_tracker)
df.to_csv('../../mTBI_scRNA_seq/DEGs/up&down_degs_across137days.csv', index=False)
