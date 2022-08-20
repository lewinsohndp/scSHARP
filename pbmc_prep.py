import pandas as pd
import numpy as np
from scipy.io import mmread
import random

"""Script to read in and combine all pbmcs, then split"""
output_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_10x_new/"
data_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_10x/"
cell_types = ["b_cells", "cd14_monocytes", "cd34","cd4_t_helper", "cd56_nk", "cytotoxic_t", "memory_t", "naive_cytotoxic", "naive_t", "regulatory_t"]

random.seed(8)
random_suffixes=random.sample(range(10,100),len(cell_types))
final_df = None
labels = []
for i, cell_type in enumerate(cell_types):
    temp_data = mmread(data_path + cell_type + "_filtered_gene_bc_matrices/hg19/matrix.mtx").toarray()
    temp_genes = pd.read_csv(data_path + cell_type + "_filtered_gene_bc_matrices/hg19/genes.tsv", sep="/t", header=None).to_numpy().flatten()
    temp_cells = pd.read_csv(data_path + cell_type + "_filtered_gene_bc_matrices/hg19/barcodes.tsv", sep="/t", header=None).to_numpy().flatten()
        
    if i == 0:
        final_df = pd.DataFrame(temp_data)
        final_df.index = temp_genes
        final_df.columns = temp_cells
        final_df = final_df.add_suffix("-" + str(random_suffixes[i]))
    else:
        temp_df = pd.DataFrame(temp_data)
        temp_df.index = temp_genes
        temp_df.columns = temp_cells
        temp_df = temp_df.add_suffix("-" + str(random_suffixes[i]))
        final_df = final_df.merge(temp_df, how='inner', left_index=True, right_index=True)
    
    labels += [cell_type]*len(temp_cells)
    print(temp_data.shape)
print(final_df.shape)
final_df = final_df.T
print(final_df.shape)

final_df['labels'] = labels

shuffled_df = final_df.sample(frac=1, random_state=8)
print(shuffled_df.shape)

split_one = shuffled_df.sample(frac=0.10, random_state=9)
split_two = shuffled_df.drop(split_one.index)

print(split_one.shape)
print(split_two.shape)

split_one_ref = split_one.sample(frac=0.5, random_state=10)
split_one_query = split_one.drop(split_one_ref.index)

print(split_one_ref.shape)
print(split_one_query.shape)

split_one_ref_labels = split_one_ref['labels']
split_one_query_labels = split_one_query['labels']
split_two_labels = split_two['labels']

split_one_ref = split_one_ref.drop("labels", axis=1)
split_one_query = split_one_query.drop("labels", axis=1)
split_two = split_two.drop("labels", axis=1)

split_one_ref.to_csv(output_path + "smaller_split_ref.csv")
split_one_query.to_csv(output_path + "smaller_split_query.csv")
split_one_ref_labels.to_csv(output_path + "smaller_split_labels_ref.csv")
split_one_query_labels.to_csv(output_path + "smaller_split_labels_query.csv")

split_two.to_csv(output_path + "bigger_split.csv")
split_two_labels.to_csv(output_path + "bigger_split_labels.csv")

