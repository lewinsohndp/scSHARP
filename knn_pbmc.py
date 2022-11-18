import torch
import utilities
import numpy as np
import pandas as pd
from sc_sharp import scSHARP

data_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/"

counts = pd.read_csv(data_path + "counts.csv", index_col = 0)
all_labels = pd.read_csv(data_path + "preds.csv", index_col=0)
X, keep_cells,keep_genes,pca_obj = utilities.preprocess(np.array(counts), scale=False, run_pca=True)

_,marker_names = utilities.read_marker_file(data_path + "markers.txt")
all_labels_factored = utilities.factorize_df(all_labels, marker_names)
encoded_labels = utilities.encode_predictions(all_labels_factored)
confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = .51)
confident_labels = confident_labels
train_nodes = np.where(confident_labels != -1)[0]
test_nodes = np.where(confident_labels == -1)[0]

real_y = pd.read_csv(data_path + "labels_cd4-8.csv", index_col= 0).iloc[:,0]
real_y = pd.factorize(real_y, sort=True)[0]
real_y = real_y[keep_cells]

data_name = "pbmc"
neighbors = [10,50,100, 200, 300, 400, 500, 600, 700]
accuracies = []
# one epoch
preds = utilities.knn_consensus_batch(X, confident_labels, 200, converge=True, one_epoch=False, batch_size=1000, keep_conf=True)
print(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), train_nodes, test_nodes))

