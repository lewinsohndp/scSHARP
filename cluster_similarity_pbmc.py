import utilities
import numpy as np
import pandas as pd
from itertools import combinations

def compare_cell_clusters(data_folder, counts, meta, meta_col):
    # read in dataset
    data_path = data_folder + counts
    X = pd.read_csv(data_path, index_col=0)
    X, keep_cells,_,_ = utilities.preprocess(np.array(X), scale=False)
    meta_path = data_folder + meta
    metadata = pd.read_csv(meta_path, index_col=0)
    if type(meta_col)==int:
        real_y, key = pd.factorize(metadata.iloc[:,meta_col], sort=True)
    else: real_y, key = pd.factorize(metadata[meta_col], sort=True)
    real_y = real_y[keep_cells]
    
    for combo in combinations(np.unique(real_y),2):
        disimilarity = utilities.average_linkage(X[real_y == combo[0]], X[real_y == combo[1]])
        print(str(key[combo[0]]) + " and " + str(key[combo[1]]) + " = " + str(disimilarity))
    
    print(key)
    return X, real_y

X, real_y = compare_cell_clusters("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/", 'counts.csv', "labels_cd4-8.csv", 0)
