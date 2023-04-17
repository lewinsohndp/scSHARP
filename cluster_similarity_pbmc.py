import utilities
import numpy as np
import pandas as pd
from itertools import combinations, combinations_with_replacement
from sklearn.metrics import pairwise_distances

def disimilarity_new(X, group_1, group_2):
    group_1 = np.array(group_1)
    group_2 = np.array(group_2)
    dist = pairwise_distances(X[group_1,:], X[group_2,:])
    total = np.sum(dist)
    return total / (group_1.sum() * group_2.sum())
    
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
    
    print(len(real_y))
    print(X.shape)
    types = len(np.unique(real_y))
    ret = np.zeros(shape=(types, types))
    for combo in combinations_with_replacement(np.unique(real_y),2):
        #disimilarity = utilities.average_linkage(X[real_y == combo[0]], X[real_y == combo[1]])
        disimilarity = disimilarity_new(X, real_y == combo[0], real_y == combo[1])
        ret[combo[0], combo[1]] = disimilarity
        print(str(key[combo[0]]) + " and " + str(key[combo[1]]) + " = " + str(disimilarity))
        
    upper = np.triu(ret)
    ret = upper + upper.T - np.diag(upper.diagonal())
    
    print(key)
    return ret

mat = compare_cell_clusters("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/", 'counts.csv', "labels_cd4-8.csv", 0)

print(mat)
