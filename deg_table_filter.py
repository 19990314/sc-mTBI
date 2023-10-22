import os
import pandas as pd
import logging

# Set up the logger
logging.basicConfig(filename='mylog.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.critical('\n======================================New Pipeline Started======================================')




# global variabless

# output_directory
output_dir = "../../mTBI_scRNA_seq/Oct20/"

# thresholds
p_adj = 5e-12 #
pv = 0.001
fc = 0.25

# statistics
temporal_deg_tracker = {
    'studies': ["Sham v.s. Day1", "Sham v.s. Day3", "Sham v.s. Day7", "Sham v.s. Combined"],
    'upregulated': [],
    'up_count': [],
    'downregulated': [],
    'down_count': []
}

celltype_deg_tracker = {
    'celltypes': ["aEC", "vEC", "capEC"],
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
    'celltypes': ["aEC", "vEC", "capEC", "allEC"],
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

def read_deg_tables(filename_prefix):
    # Directory containing the CSV files
    directory = output_dir + "original_DEGs/"

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


def filter_degs(dataframes, pval, fc, pval_adj):
    for celltype_key in dataframes:
        size_tracker['celltype'].append(celltype_key)

        # before filtering
        size_tracker['before filtering'].append(len(dataframes[celltype_key]))
        # check Tmem252
        if "Tmem252" in dataframes[celltype_key].index:
            size_tracker['Tmem252 check1'].append("Yes")
        else:
            size_tracker['Tmem252 check1'].append("No")

        # filter by p val
        dataframes[celltype_key] = dataframes[celltype_key][dataframes[celltype_key]['p_val'] <= pval]
        size_tracker['after p-val'].append(len(dataframes[celltype_key]))
        # check Tmem252
        if "Tmem252" in dataframes[celltype_key].index:
            size_tracker['Tmem252 check2'].append("Yes")
        else:
            size_tracker['Tmem252 check2'].append("No")

        # filter by fc
        dataframes[celltype_key] = dataframes[celltype_key][abs(dataframes[celltype_key]['avg_log2FC']) >= fc]
        size_tracker['after fc'].append(len(dataframes[celltype_key]))
        #check Tmem252
        if "Tmem252" in dataframes[celltype_key].index:
            size_tracker['Tmem252 check3'].append("Yes")
        else:
            size_tracker['Tmem252 check3'].append("No")


        # filter by adj p-val
        dataframes[celltype_key] = dataframes[celltype_key][dataframes[celltype_key]['p_val_adj'] <= pval_adj]
        size_tracker['after p-adj'].append(len(dataframes[celltype_key]))
        #check Tmem252
        if "Tmem252" in dataframes[celltype_key].index:
            size_tracker['Tmem252 check4'].append("Yes")
        else:
            size_tracker['Tmem252 check4'].append("No")
    return dataframes


def divide_up_down_regulation(dataframes, day_index, pct):
    logging.info('Up- and down- regulated DEGs for each cell type were saved to: DEGs_filtered_updown_sorted/')  # *log

    for celltype_key in dataframes:
        temp_dic = {}

        # extract up-regulated, plus sort by fc
        temp_dic["up"] = dataframes[celltype_key][dataframes[celltype_key]['avg_log2FC'] > 0].sort_values(by='avg_log2FC', ascending=False)

        # filter by pct
        temp_dic["up"] = temp_dic["up"][temp_dic["up"]['pct.1'] >= pct]

        # output
        temp_dic["up"].to_csv(output_dir + 'DEGs_filtered_updown_sorted/DEG_' + day_index + "_" + celltype_key + "_up_sorted.csv") # output
        logging.critical('Output ---> DEGs_filtered_updown_sorted/DEG_' + day_index + "_" + celltype_key + "_up_sorted.csv")  # *log

        # extract down-regulated, plus sort by fc
        temp_dic["down"] = dataframes[celltype_key][dataframes[celltype_key]['avg_log2FC'] < 0].sort_values(by='avg_log2FC', ascending=True)

        # filter by pct
        temp_dic["down"] = temp_dic["down"][temp_dic["down"]['pct.2'] >= pct]

        # output
        temp_dic["down"].to_csv(output_dir + 'DEGs_filtered_updown_sorted/DEG_' + day_index + "_" + celltype_key + "_down_sorted.csv") # output
        logging.critical('Output ---> DEGs_filtered_updown_sorted/DEG_' + day_index + "_" + celltype_key + "_down_sorted.csv")  # *log

        # replace unsplitted df
        dataframes[celltype_key] = temp_dic

        # statistics
        size_tracker['after pct (UP/DOWN)'].append(str(len(dataframes[celltype_key]["up"])) + "/" + str(len(dataframes[celltype_key]["down"])))

    return dataframes

def celltype_specific_DEGs(dataframes, day_index, up, down):
    logging.info('Up- and down- regulated DEGs (celltype specific) were saved to: DEGs_celltype_specific/')  # *log

    for celltype_key in dataframes:
        # up
        dataframes[celltype_key]["up"].drop(set(up).intersection(dataframes[celltype_key]["up"].index)).to_csv(
            output_dir + 'DEGs_celltype_specific/DEG_' + day_index + "_" + celltype_key + "_up_sorted.csv")  # output
        logging.critical(
            'Output ---> DEGs_celltype_specific/DEG_' + day_index + "_" + celltype_key + "_up_sorted.csv")


        # down
        dataframes[celltype_key]["down"].drop(set(down).intersection(dataframes[celltype_key]["down"].index)).to_csv(
            output_dir + 'DEGs_celltype_specific/DEG_' + day_index + "_" + celltype_key + "_down_sorted.csv")  # output
        logging.critical(
            'Output ---> DEGs_celltype_specific/DEG_' + day_index + "_" + celltype_key + "_down_sorted.csv")



    return dataframes


def find_common_degs_for_a_day(dataframes):
    common_degs_up = {}
    common_degs_down = {}
    unique_up = ()
    unique_down = ()

    # dataframes: a list of deg dfs: aEC, vEC, capEC, allEC
    for df_key in dataframes:
        # up DEGs
        for gene, info in dataframes[df_key]["up"].iterrows():
            if gene not in common_degs_up.keys():
                common_degs_up[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_up[gene] += abs(info['avg_log2FC'])

        # down DEGs
        for gene, info in dataframes[df_key]["down"].iterrows():
            if gene not in common_degs_down.keys():
                common_degs_down[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_down[gene] += abs(info['avg_log2FC'])

        if len(unique_up) == 0:
            unique_up = set(dataframes[df_key]["up"].index)
        else:
            unique_up = unique_up.intersection(dataframes[df_key]["up"].index)

        if len(unique_down) == 0:
            unique_down = set(dataframes[df_key]["down"].index)
        else:
            unique_down = unique_down.intersection(dataframes[df_key]["down"].index)


    # Sort the DEGs
    sorted_common_degs_up = dict(sorted(common_degs_up.items(), key=lambda kv: kv[1], reverse=True))
    sorted_common_degs_down = dict(sorted(common_degs_down.items(), key=lambda kv: kv[1], reverse=True))

    common_degs_up = set(sorted_common_degs_up.keys())
    common_degs_down = set(sorted_common_degs_down.keys())

    # check appear twice
    for gene in common_degs_up:
        score = 0
        for df_key in dataframes:
            if gene in dataframes[df_key]["up"].index:
                score += 1
        if score < 2:
            sorted_common_degs_up.pop(gene)

    for gene in common_degs_down:
        score = 0
        for df_key in dataframes:
            if gene in dataframes[df_key]["down"].index:
                score += 1
        if score < 2:
            sorted_common_degs_down.pop(gene)


    return list(sorted_common_degs_up.keys()), list(sorted_common_degs_down.keys()), list(unique_up), list(unique_down)


def find_common_degs_for_a_celltype(day1_df, day3_df, day7_df):
    up = []
    down = []
    up_ct = []
    down_ct = []

    # dataframes: a list of deg dfs: aEC, vEC, capEC, allEC
    for celltype_key in ["aEC", "vEC", "capEC"]:
        common_degs_up = {}
        common_degs_down = {}

        for gene, info in day1_df[celltype_key]["up"].iterrows():
            if gene not in common_degs_up.keys():
                common_degs_up[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_up[gene] += abs(info['avg_log2FC'])

        for gene, info in day1_df[celltype_key]["down"].iterrows():
            if gene not in common_degs_down.keys():
                common_degs_down[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_down[gene] += abs(info['avg_log2FC'])

        for gene, info in day3_df[celltype_key]["up"].iterrows():
            if gene not in common_degs_up.keys():
                common_degs_up[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_up[gene] += abs(info['avg_log2FC'])

        for gene, info in day3_df[celltype_key]["down"].iterrows():
            if gene not in common_degs_down.keys():
                common_degs_down[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_down[gene] += abs(info['avg_log2FC'])

        for gene, info in day7_df[celltype_key]["up"].iterrows():
            if gene not in common_degs_up.keys():
                common_degs_up[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_up[gene] += abs(info['avg_log2FC'])

        for gene, info in day7_df[celltype_key]["down"].iterrows():
            if gene not in common_degs_down.keys():
                common_degs_down[gene] = abs(info['avg_log2FC'])
            else:
                common_degs_down[gene] += abs(info['avg_log2FC'])

        # Sort the DEGs
        sorted_common_degs_up = dict(sorted(common_degs_up.items(), key=lambda kv: kv[1], reverse=True))
        sorted_common_degs_down = dict(sorted(common_degs_down.items(), key=lambda kv: kv[1], reverse=True))

        #
        up.append(list(sorted_common_degs_up.keys()))
        up_ct.append(len(sorted_common_degs_up.keys()))
        down.append(list(sorted_common_degs_down.keys()))
        down_ct.append(len(sorted_common_degs_down.keys()))

    # Find intersections
    intersection_all_celltypes_up = list(set(up[0]).intersection(up[1], up[2]))
    intersection_all_celltypes_down = list(set(down[0]).intersection(down[1], down[2]))

    up_celltype_specific = up
    down_celltype_specific = down
    up_celltype_specific_ct = []
    down_celltype_specific_ct = []

    # Remove elements that are in intersections
    for item in intersection_all_celltypes_up:
        for celltype_deg in up_celltype_specific:
            if item in celltype_deg:
                celltype_deg.remove(item)
    for i in up_celltype_specific:
        up_celltype_specific_ct.append(len(i))


    for item in intersection_all_celltypes_down:
        for celltype_deg in down_celltype_specific:
            if item in celltype_deg:
                celltype_deg.remove(item)
    for i in down_celltype_specific:
        down_celltype_specific_ct.append(len(i))

    return up, up_ct, down, down_ct, up_celltype_specific, up_celltype_specific_ct, down_celltype_specific, down_celltype_specific_ct


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



# create dir
# DEGs_filtered_updown_sorted/   filtration_statistics/   common_DEGs/
prep_workspace(output_dir + "DEGs_filtered_updown_sorted/")
prep_workspace(output_dir + "filtration_statistics/")
prep_workspace(output_dir + "Fig1g_common_DEGs/")
prep_workspace(output_dir + "DEGs_celltype_specific/")
prep_workspace(output_dir + "DEGs_day_specific/")

# *log
logging.info('Workspace Created: DEGs_filtered_updown_sorted/ filtration_statistics/ Fig1g_common_DEGs/')
logging.info('*DEGs_filtered_updown_sorted/: DEG tables for each single day and each cell types')
logging.info('*filtration_statistics/: the DEG counts and Tmem252 flag after each filtration step (each single day, each cell types)')
logging.info('*Fig1g_common_DEGs/: common DEGs and counts (each single day, each cell types)')
logging.info('*DEGs_celltype_specific/: DEGs for each cell type')

## ---------day1---------
logging.critical('********************Day 1************************')# *log

# DEG size tracker
size_tracker = {
        'celltype': [],
        'before filtering': [],
        'Tmem252 check1': [],
        'after p-val': [],
        'Tmem252 check2': [],
        'after fc': [],
        'Tmem252 check3': [],
        'after p-adj': [],
        'Tmem252 check4': [],
        'after pct (UP/DOWN)': []
}

# read
logging.info('Day 1: Read 4 DEG files (aEC, vEC, capEC, allEC)') # *log
degs_for_day1 = read_deg_tables("DEG_day1") # INPUT modify here
logging.info('---- DONE ----') # *log

# filter
logging.info('DEG filtration (p-value: 0.05, fc: 0.2, pval-adj: 0.001)')# *log
degs_for_day1 = filter_degs(degs_for_day1, pv, fc,  p_adj)
logging.info('---- DONE ----') # *log


# up & down
logging.info('DEG UP or Down (pct: 0.1)')# *log
degs_for_day1 = divide_up_down_regulation(degs_for_day1, "day1", 0.1)
logging.info('---- DONE ----') # *log


# diff degs
logging.info('Find common DEGs for Day1 up- and down-regulated DEGs from its a,v,cap,allEC groups')# *log
up, down, up_unique, down_unique = find_common_degs_for_a_day(degs_for_day1)
logging.info('---- DONE ----') # *log

# cell-type specific DEG tables
logging.info('Find cell-type specific DEGs for Day1 up- and down-regulated DEGs from its a,v,cap groups')# *log
celltype_specific_DEGs(degs_for_day1, "day1", up, down)
logging.info('---- DONE ----') # *log

# record day1 DEGs
temporal_deg_tracker['upregulated'].append(up)
temporal_deg_tracker['up_count'].append(len(up))
temporal_deg_tracker['downregulated'].append(down)
temporal_deg_tracker['down_count'].append(len(down))

# output
df = pd.DataFrame(size_tracker)
df.to_csv(output_dir + 'filtration_statistics/DEG_statistics_day1.csv', index=False) # OUTPUT modify here
logging.critical('OUTPUT ---> filtration_statistics/DEG_statistics_day1.csv')# *log






## ---------day3---------
logging.critical('********************Day 3************************')# *log

# DEG size tracker
size_tracker = {
    'celltype': [],
    'before filtering': [],
    'Tmem252 check1': [],
    'after p-val': [],
    'Tmem252 check2': [],
    'after fc': [],
    'Tmem252 check3': [],
    'after p-adj': [],
    'Tmem252 check4': [],
    'after pct (UP/DOWN)': []

}

# read
logging.info('Day 3: Read 4 DEG files (aEC, vEC, capEC, allEC)') # *log
degs_for_day3 = read_deg_tables("DEG_day3")
logging.info('---- DONE ----') # *log

# filter
logging.info('DEG filtration (p-value: 0.05, fc: 0.2, pval-adj: 0.001)')# *log
degs_for_day3 = filter_degs(degs_for_day3, pv, fc,  p_adj)
logging.info('---- DONE ----') # *log


# up & down
logging.info('DEG UP or Down (pct: 0.1)')# *log
degs_for_day3 = divide_up_down_regulation(degs_for_day3, "day3", 0.1)
logging.info('---- DONE ----')# *log


# diff degs
logging.info('Find common DEGs for Day3 up- and down-regulated DEGs from its a,v,cap,allEC groups')# *log
up, down, up_unique, down_unique = find_common_degs_for_a_day(degs_for_day3)
logging.info('---- DONE ----')# *log


# cell-type specific DEG tables
logging.info('Find cell-type specific DEGs for Day1 up- and down-regulated DEGs from its a,v,cap groups')# *log
celltype_specific_DEGs(degs_for_day3, "day3", up, down)
logging.info('---- DONE ----') # *log


# record day1 DEGs
temporal_deg_tracker['upregulated'].append(up)
temporal_deg_tracker['up_count'].append(len(up))
temporal_deg_tracker['downregulated'].append(down)
temporal_deg_tracker['down_count'].append(len(down))

# output
df = pd.DataFrame(size_tracker)
df.to_csv(output_dir + 'filtration_statistics/DEG_statistics_day3.csv', index=False) # OUTPUT modify here
logging.critical('OUTPUT ---> filtration_statistics/DEG_statistics_day3.csv')# *log




## ---------day7---------
logging.critical('********************Day 7************************')# *log

# DEG size tracker
size_tracker = {
    'celltype': [],
    'before filtering': [],
    'Tmem252 check1': [],
    'after p-val': [],
    'Tmem252 check2': [],
    'after fc': [],
    'Tmem252 check3': [],
    'after p-adj': [],
    'Tmem252 check4': [],
    'after pct (UP/DOWN)': []

}

# read
logging.info('Day 3: Read 4 DEG files (aEC, vEC, capEC, allEC)') # *log
degs_for_day7 = read_deg_tables("DEG_day7")
logging.info('---- DONE ----') # *log

# filter
logging.info('DEG filtration (p-value: 0.05, fc: 0.2, pval-adj: 0.001)')# *log
degs_for_day7 = filter_degs(degs_for_day7, pv, fc,  p_adj)
logging.info('---- DONE ----') # *log


# up & down
logging.info('DEG UP or Down (pct: 0.1)')# *log
degs_for_day7 = divide_up_down_regulation(degs_for_day7, "day7", 0.1)
logging.info('---- DONE ----')# *log


# diff degs
logging.info('Find common DEGs for Day7 up- and down-regulated DEGs from its a,v,cap,allEC groups')# *log
up, down, up_unique, down_unique = find_common_degs_for_a_day(degs_for_day7)
logging.info('---- DONE ----')# *log


# cell-type specific DEG tables
logging.info('Find cell-type specific DEGs for Day1 up- and down-regulated DEGs from its a,v,cap groups')# *log
celltype_specific_DEGs(degs_for_day7, "day7", up, down)
logging.info('---- DONE ----') # *log


# record day1 DEGs
temporal_deg_tracker['upregulated'].append(up)
temporal_deg_tracker['up_count'].append(len(up))
temporal_deg_tracker['downregulated'].append(down)
temporal_deg_tracker['down_count'].append(len(down))


# output
df = pd.DataFrame(size_tracker)
df.to_csv(output_dir + 'filtration_statistics/DEG_statistics_day7.csv', index=False) # OUTPUT modify here
logging.critical('OUTPUT ---> filtration_statistics/DEG_statistics_day7.csv')# *log


## ---------all days combined---------
logging.critical('********************All days************************')# *log

# DEG size tracker
size_tracker = {
    'celltype': [],
    'before filtering': [],
    'Tmem252 check1': [],
    'after p-val': [],
    'Tmem252 check2': [],
    'after fc': [],
    'Tmem252 check3': [],
    'after p-adj': [],
    'Tmem252 check4': [],
    'after pct (UP/DOWN)': []
}

# read
logging.info('Day 3: Read 4 DEG files (aEC, vEC, capEC, allEC)') # *log
degs_for_alldays = read_deg_tables("DEG_combined")
logging.info('---- DONE ----') # *log

# filter
logging.info('DEG filtration (p-value: 0.05, fc: 0.2, pval-adj: 0.001)')# *log
degs_for_alldays = filter_degs(degs_for_alldays, pv, fc,  p_adj)
logging.info('---- DONE ----') # *log


# up & down
logging.info('DEG UP or Down (pct: 0.1)')# *log
degs_for_alldays = divide_up_down_regulation(degs_for_alldays, "combined", 0.1)
logging.info('---- DONE ----')# *log


# diff degs
logging.info('Find common DEGs for Day1 up- and down-regulated DEGs from its a,v,cap,allEC groups')# *log
up, down, up_unique, down_unique = find_common_degs_for_a_day(degs_for_alldays)
logging.info('---- DONE ----')# *log


# cell-type specific DEG tables
logging.info('Find cell-type specific DEGs for Day1 up- and down-regulated DEGs from its a,v,cap groups')# *log
celltype_specific_DEGs(degs_for_alldays, "alldays", up, down)
logging.info('---- DONE ----') # *log


# record day1 DEGs
temporal_deg_tracker['upregulated'].append(up)
temporal_deg_tracker['up_count'].append(len(up))
temporal_deg_tracker['downregulated'].append(down)
temporal_deg_tracker['down_count'].append(len(down))

# output
df = pd.DataFrame(size_tracker)
df.to_csv(output_dir + 'filtration_statistics/DEG_statistics_alldays.csv', index=False) # OUTPUT modify here
logging.critical('OUTPUT ---> filtration_statistics/DEG_statistics_alldays.csv')# *log






## ---------find common deg for each cell type---------
celltype_deg_tracker["upregulated"], celltype_deg_tracker['up_count'], celltype_deg_tracker["downregulated"], celltype_deg_tracker['down_count'], celltype_deg_tracker['upregulated_celltype_specific'], celltype_deg_tracker['up_celltype_specific_count'], celltype_deg_tracker['downregulated_celltype_specific'], celltype_deg_tracker['down_celltype_specific_count'] = find_common_degs_for_a_celltype(degs_for_day1, degs_for_day3, degs_for_day7)


df = pd.DataFrame(celltype_deg_tracker)
df.to_csv(output_dir + "Fig1g_common_DEGs/up&down_degs_celltype_specific.csv", index=False) # OUTPUT modify here
logging.info('Output ---> Fig1g_common_DEGs/up&down_degs_celltype_specific.csv')  # *log

# --------- output degs for day 1,3,7 ---------
df = pd.DataFrame(temporal_deg_tracker)
df.to_csv(output_dir + "Fig1g_common_DEGs/up&down_degs_temporal_specific.csv", index=False) # OUTPUT modify here
logging.info('Output ---> Fig1g_common_DEGs/up&down_degs_temporal_specific.csv')  # *log



## ---------temporal-specific DEGs---------
union_DEG = set(degs_for_day1["allEC"]["up"].index).union(degs_for_day3["allEC"]["up"].index, degs_for_day7["allEC"]["up"].index)
common_DEG_alldays = []
for i in union_DEG:
    score = 0
    if i in degs_for_day1["allEC"]["up"].index:
        score += 1
    if i in degs_for_day3["allEC"]["up"].index:
        score += 1
    if i in degs_for_day7["allEC"]["up"].index:
        score += 1

    if score>=2:
        common_DEG_alldays.append(i)

degs_for_day1["allEC"]["up"].drop(set(common_DEG_alldays).intersection(degs_for_day1["allEC"]["up"].index)).to_csv(output_dir + 'DEGs_day_specific/DEG_1_allEC_up_sorted.csv')
degs_for_day3["allEC"]["up"].drop(set(common_DEG_alldays).intersection(degs_for_day3["allEC"]["up"].index)).to_csv(output_dir + 'DEGs_day_specific/DEG_3_allEC_up_sorted.csv')
degs_for_day7["allEC"]["up"].drop(set(common_DEG_alldays).intersection(degs_for_day7["allEC"]["up"].index)).to_csv(output_dir + 'DEGs_day_specific/DEG_7_allEC_up_sorted.csv')


union_DEG = set(degs_for_day1["allEC"]["down"].index).union(degs_for_day3["allEC"]["down"].index, degs_for_day7["allEC"]["down"].index)
common_DEG_alldays = []
for i in union_DEG:
    score = 0
    if i in degs_for_day1["allEC"]["down"].index:
        score += 1
    if i in degs_for_day3["allEC"]["down"].index:
        score += 1
    if i in degs_for_day7["allEC"]["down"].index:
        score += 1

    if score >= 2:
        common_DEG_alldays.append(i)

degs_for_day1["allEC"]["down"].drop(set(common_DEG_alldays).intersection(degs_for_day1["allEC"]["down"].index)).to_csv(output_dir + 'DEGs_day_specific/DEG_1_allEC_down_sorted.csv')
degs_for_day3["allEC"]["down"].drop(set(common_DEG_alldays).intersection(degs_for_day3["allEC"]["down"].index)).to_csv(output_dir + 'DEGs_day_specific/DEG_3_allEC_down_sorted.csv')
degs_for_day7["allEC"]["down"].drop(set(common_DEG_alldays).intersection(degs_for_day7["allEC"]["down"].index)).to_csv(output_dir + 'DEGs_day_specific/DEG_7_allEC_down_sorted.csv')



