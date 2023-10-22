import pandas as pd
import os


# global variables

# output_directory
output_dir = "../../mTBI_scRNA_seq/Oct17/"
pathway_dir = output_dir + "pathway/"


def read_pathway_tables(filename_prefix):

    # Get a list of CSV files with the specified prefix
    excel_files = [file for file in os.listdir(pathway_dir) if (file.startswith(filename_prefix)) & ("sort" not in file)]

    # Create an empty dictionary to store dataframes
    dataframes = {}
    dataframes = []
    merged = None

    # Loop through the CSV files and read them into dataframes
    for excel_file in excel_files:
        file_path = os.path.join(pathway_dir, excel_file)
        # Use the file name (without extension) as the dataframe key
        #dataframe_key = get_EC_celltype(os.path.splitext(excel_file)[0])

        # Read the CSV file into a dataframe
        #dataframes[dataframe_key] = pd.read_csv(file_path, header=0, names = ['gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj'], dtype='object')
        #dataframes[dataframe_key].set_index('gene', inplace=True)
        #dataframes[dataframe_key] = dataframes[dataframe_key].apply(pd.to_numeric, errors='coerce')

        current_df = pd.read_excel(file_path, header=True)

        if merged == None:
            merged = current_df
        else:
            merged = pd.merge(merged, current_df, on='key')


    return merged

merged = read_pathway_tables("pathway")
print(merged)