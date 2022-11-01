import torch
import utilities
import numpy as np
import pandas as pd
from sc_sharp import scSHARP

data_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/"

counts = pd.read_csv(data_path + "counts.csv", index_col = 0)
all_labels = pd.read_csv(data_path + "preds.csv", index_col=0)
X, keep_cells,keep_genes,pca_obj = utilities.preprocess(np.array(counts), scale=False, run_pca=False)

_,marker_names = utilities.read_marker_file(data_path + "markers.txt")
all_labels_factored = utilities.factorize_df(all_labels, marker_names)
encoded_labels = utilities.encode_predictions(all_labels_factored)
confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = .51)
confident_labels = confident_labels

real_y = pd.read_csv(data_path + "labels_cd4-8.csv", index_col= 0).iloc[:,0]
real_y = pd.factorize(real_y, sort=True)[0]
real_y = real_y[keep_cells]

print(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(confident_labels), range(len(real_y)), range(len(real_y)))[0:2])

# one epoch
for i in [10,50,100,200]:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=False, one_epoch=True, batch_size=1000)
    print(i)
    print(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), range(len(real_y)), range(len(real_y)))[0:2])


# until converge
for i in [10,50,100,200]:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=True, one_epoch=False, batch_size=1000)
    print(i)
    print(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), range(len(real_y)), range(len(real_y)))[0:2])

# until all labelled
for i in [10,50,100,200]:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=False, one_epoch=False, batch_size=1000)
    print(i)
    print(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), range(len(real_y)), range(len(real_y)))[0:2])


