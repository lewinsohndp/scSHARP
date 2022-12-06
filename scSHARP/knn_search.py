import sys
import numpy as np
import pandas as pd
import json
import torch
from . import utilities
from .gcn_model import GCNModel
import os
import statistics
import random

"""script that tests knn method for given dataset"""

# for extended abstract, keep_conf = False used
keep_confident = True
data_folder = sys.argv[1]
marker_path = data_folder + "markers.txt"
if data_folder == "/home/groups/ConradLab/daniel/sharp_data/jung/":
    tools = ["sctype","scsorter","scina"]
else:
    tools = ["sctype","scsorter","scina", "singler", "scpred"]
votes = .51

print(data_folder)
all_labels = pd.read_csv(data_folder + "preds.csv", index_col=0)
if all_labels.shape[1] != len(tools): raise Exception("wrong amount of tools in file")

data_path = data_folder + "query_counts.csv"
#data_path = data_folder + "counts.csv"

# read in dataset
X = pd.read_csv(data_path, index_col=0)
X, keep_cells,_,_ = utilities.preprocess(np.array(X), scale=False, run_pca=True)

_,marker_names = utilities.read_marker_file(marker_path)

all_labels_factored = utilities.factorize_df(all_labels, marker_names)
encoded_labels = utilities.encode_predictions(all_labels_factored)

confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = votes)
train_nodes = np.where(confident_labels != -1)[0]
test_nodes = np.where(confident_labels == -1)[0]

#create validation set
random.seed(8)
np.random.seed(8)
#validation_nodes = random.sample(train_nodes, int(len(train_nodes)*.2))
validation_nodes = np.random.choice(train_nodes, size=int(len(train_nodes)*.2), replace=False)
unmasked_confident = np.copy(confident_labels)
confident_labels[validation_nodes] = -1

print("One epoch")
#one epoch
for i in [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900]:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=False, one_epoch=True, batch_size=1000, keep_conf = keep_confident)
    print(i)
    print(utilities.validation_metrics(torch.tensor(unmasked_confident), torch.tensor(preds), validation_nodes, range(len(unmasked_confident)))[0:4])

print("Until all labelled")
# until all labelled
for i in [10,50,100, 200, 300, 400, 500, 600, 700, 800, 900]:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=False, one_epoch=False, batch_size=1000, keep_conf = keep_confident)
    print(i)
    print(utilities.validation_metrics(torch.tensor(unmasked_confident), torch.tensor(preds), validation_nodes, range(len(unmasked_confident)))[0:4])

print("Until Converge")
# until converge batch
for i in [10,50,100,200,300, 400,500,600,700,800,900]:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, batch_size = 1000, converge=True, one_epoch=False, keep_conf = keep_confident)
    print(i)
    print(utilities.validation_metrics(torch.tensor(unmasked_confident), torch.tensor(preds), validation_nodes, range(len(unmasked_confident)))[0:4])



