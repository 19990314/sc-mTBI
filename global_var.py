import os
import pandas as pd

# global variabless

# output_directory
output_dir = "../../mTBI_scRNA_seq/Jan12/"
updown_DEG_dir = output_dir + "DEGs_filtered_updown_sorted/"

# thresholds
p_adj = 0.05 #
pv = 0.05
fc = 0.25
chaos_pct = 0.5

cell_types = ["aEC", "vEC", "capEC","acapEC","vcapEC"]

# statistics
temporal_deg_tracker = {
    'studies': ["Sham v.s. Day1", "Sham v.s. Day3", "Sham v.s. Day7", "Sham v.s. Combined"],
    'upregulated': [],
    'up_count': [],
    'downregulated': [],
    'down_count': []
}


celltype_deg_tracker = {
    'celltypes': cell_types,
    'upregulated': [],
    'up_count': [],
    'downregulated': [],
    'down_count': [],
    'upregulated_celltype_specific': [],
    'up_celltype_specific_count': [],
    'downregulated_celltype_specific': [],
    'down_celltype_specific_count': []
}

celltype_specific_deg_tracker = {
    'celltypes': cell_types.append("allEC"),
    'upregulated': [],
    'up_count': [],
    'downregulated': [],
    'down_count': []
}



def prep_workspace(dir):
    # Check if the directory exists
    if not os.path.exists(dir):
        os.mkdir(dir)


def get_EC_celltype(filename):
    for i in filename.split("_"):
        if "EC" in i:
            return i


def read_deg_tables(dir, filename_prefix):
    # Directory containing the CSV files
    directory = output_dir + dir

    # Get a list of CSV files with the specified prefix
    csv_files = [file for file in os.listdir(directory) if (file.startswith(filename_prefix)) & ("sort" not in file)]

    # Create an empty dictionary to store dataframes
    dataframes = {}

    # Loop through the CSV files and read them into dataframes
    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        # Use the file name (without extension) as the dataframe key
        dataframe_key = get_EC_celltype(os.path.splitext(csv_file)[0])
        # Read the CSV file into a dataframe
        dataframes[dataframe_key] = pd.read_csv(file_path, header=0, names = ['gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj'], dtype='object')
        dataframes[dataframe_key].set_index('gene', inplace=True)
        dataframes[dataframe_key] = dataframes[dataframe_key].apply(pd.to_numeric, errors='coerce')

    return dataframes
