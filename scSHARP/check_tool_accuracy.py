import numpy as np
import pandas as pd
import torch
from . import utilities
from .sc_sharp import scSHARP
import os

def get_class_accuracy(conf_mat, labels):
    label_count = 0
    for i in range(conf_mat.shape[0]):
        if conf_mat.shape[0] > len(labels) and i ==0: continue
        print(labels[label_count])
        print(conf_mat[i,i] / sum(conf_mat[i,:]))
        label_count +=1

# updated marker experiment
all_labels = pd.read_csv("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/preds_att_marker_test.csv", index_col=0)
scina_preds,_ = pd.factorize(all_labels.scina, sort=True)
sctype_preds,_ = pd.factorize(all_labels.sctype, sort=True)
scsorter_preds,_ = pd.factorize(all_labels.scsorter, sort=True)

real_labels = pd.read_csv("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/preds.csv", index_col=0)
real_scina_preds,_ = pd.factorize(real_labels.scina, sort=True)
real_sctype_preds,_ = pd.factorize(real_labels.sctype, sort=True)
real_scsorter_preds,_ = pd.factorize(real_labels.scsorter, sort=True)

meta_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/labels_cd4-8.csv"
metadata = pd.read_csv(meta_path, index_col=0)
real_labels, keys = pd.factorize(metadata.iloc[:,0], sort=True)

data_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/counts.csv"
tool_preds = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/preds.csv"
tool_list = ["scina", "scsorter", "sctype", "singler", "scpred"]
#tool_list = ["scina", "scsorter", "sctype"]
marker_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/markers.txt"
neighbors=2
config="configs/2_25.txt"
sharp = scSHARP(data_path, tool_preds, tool_list, marker_path, neighbors, config)

if os.path.exists("pbmc_test_trained_model_8"):
    print("loading pre trained")
    sharp.load_model("pbmc_test_trained_model_8")
    preds, train_nodes, test_nodes, keep_cells = sharp.run_prediction(training_epochs=0, thresh=0.51, batch_size=50, seed=8)


results = utilities.validation_metrics(torch.tensor(real_labels[keep_cells]), torch.tensor(scina_preds), train_nodes, test_nodes)
print('scina results')
print(results)
print(get_class_accuracy(results[1], keys))
results_original=(utilities.validation_metrics(torch.tensor(real_labels[keep_cells]), torch.tensor(real_scina_preds), train_nodes, test_nodes))
print(results_original)
print(get_class_accuracy(results_original[1], keys))

results = utilities.validation_metrics(torch.tensor(real_labels[keep_cells]), torch.tensor(sctype_preds), train_nodes, test_nodes)
print('sctype results')
print(results)
print(get_class_accuracy(results[1], keys))
results_original=(utilities.validation_metrics(torch.tensor(real_labels[keep_cells]), torch.tensor(real_sctype_preds), train_nodes, test_nodes))
print(results_original)
print(get_class_accuracy(results_original[1], keys))

results = utilities.validation_metrics(torch.tensor(real_labels[keep_cells]), torch.tensor(scsorter_preds), train_nodes, test_nodes)
print('scsorter results')
print(results)
print(get_class_accuracy(results[1], keys))
results_original = (utilities.validation_metrics(torch.tensor(real_labels[keep_cells]), torch.tensor(real_scsorter_preds), train_nodes, test_nodes))
print(results_original)
print(get_class_accuracy(results_original[1], keys))

