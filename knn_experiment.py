import torch
import utilities
import numpy as np
import pandas as pd
from sc_sharp import scSHARP
import matplotlib.pyplot as plt

data_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/"

counts = pd.read_csv(data_path + "counts.csv", index_col = 0)
all_labels = pd.read_csv(data_path + "preds.csv", index_col=0)
X, keep_cells,keep_genes,pca_obj = utilities.preprocess(np.array(counts), scale=False, run_pca=True)

_,marker_names = utilities.read_marker_file(data_path + "markers.txt")
all_labels_factored = utilities.factorize_df(all_labels, marker_names)
encoded_labels = utilities.encode_predictions(all_labels_factored)
confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = .51)
confident_labels = confident_labels

real_y = pd.read_csv(data_path + "labels_cd4-8.csv", index_col= 0).iloc[:,0]
real_y = pd.factorize(real_y, sort=True)[0]
real_y = real_y[keep_cells]

print(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(confident_labels), range(len(real_y)), range(len(real_y)))[0:2])
data_name = "pbmc"
neighbors = [10,50,100, 200, 300, 400, 500, 600, 700, 800,900]
accuracies = []
# one epoch
for i in neighbors:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=False, one_epoch=True, batch_size=1000)
    print(i)
    accuracies.append(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), range(len(real_y)), range(len(real_y)))[0])

plt.scatter(neighbors, accuracies)
plt.xlabel("Neighbors")
plt.ylabel("Accuracy")
plt.title("One Epoch")
plt.savefig('figures/' + data_name + '_knn_plot_one.pdf')
accuracies = []
plt.clf()

# until converge
for i in neighbors:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=True, one_epoch=False, batch_size=1000)
    print(i)
    accuracies.append(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), range(len(real_y)), range(len(real_y)))[0])

plt.scatter(neighbors, accuracies)
plt.xlabel("Neighbors")
plt.ylabel("Accuracy")
plt.title("Iterate until convergence")
plt.savefig('figures/' + data_name + '_knn_plot_converge.pdf') 
accuracies = []
plt.clf()

# until all labelled
for i in neighbors:
    preds = utilities.knn_consensus_batch(X, confident_labels, i, converge=False, one_epoch=False, batch_size=1000)
    print(i)
    accuracies.append(utilities.validation_metrics(torch.tensor(real_y), torch.tensor(preds), range(len(real_y)), range(len(real_y)))[0])

plt.scatter(neighbors, accuracies)
plt.xlabel("Neighbors")
plt.ylabel("Accuracy")
plt.title("Iterate until all labeled")
plt.savefig('figures/' + data_name + '_knn_plot_all.pdf')
plt.clf()
accuracies = []
