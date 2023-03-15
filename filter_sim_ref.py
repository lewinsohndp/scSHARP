import utilities
import pandas as pd
import numpy as np

for i in [0.7, 0.8]:
    labels = pd.read_csv("/home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_" + str(i) + "_imbalance_031223/ref_labels.csv", index_col = 0)
    counts = pd.read_csv("/home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_" + str(i) +"_imbalance_031223/ref_counts.csv", index_col = 0)
    print(counts.shape)
    new_counts, row_filter,col_filter,_ = utilities.preprocess(np.array(counts), normalize=False, scale=False, run_pca=False)

    new_counts = pd.DataFrame(new_counts)
    new_counts.index = counts.index[row_filter]
    new_counts.columns = counts.columns[col_filter]
    labels = labels.iloc[row_filter,:]

    print(new_counts.shape)
    print(labels.shape)
    labels.to_csv("/home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_" + str(i) + "_imbalance_031223/ref_labels_filtered.csv")

    new_counts.to_csv("/home/groups/ConradLab/daniel/sharp_data/sharp_sims_test/splat_" + str(i) + "_imbalance_031223/ref_counts_filtered.csv")

