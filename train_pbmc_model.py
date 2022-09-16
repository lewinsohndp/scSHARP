import utilities
from sc_sharp import scSHARP
import numpy as np
import pandas as pd
import torch
import os

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
else:
    print("training new")
    preds, train_nodes, test_nodes, keep_cells = sharp.run_prediction(training_epochs=150, thresh=0.51, batch_size=50, seed=8)
    sharp.save_model("pbmc_test_trained_model_8")

meta_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/labels_cd4-8.csv"
metadata = pd.read_csv(meta_path, index_col=0)
real_y = pd.factorize(metadata.iloc[:,0], sort=True)[0]
real_y = real_y[keep_cells]
print(real_y.shape)
print(np.count_nonzero(real_y == 0))
print(np.count_nonzero(real_y == 1))
print(np.count_nonzero(real_y == 2))
print(np.count_nonzero(real_y == 3))
print(np.count_nonzero(real_y == 4))
print(np.count_nonzero(keep_cells))
print(len(test_nodes))
print(len(train_nodes))
acc, conf_mat, _,_,_,_ = utilities.validation_metrics(torch.tensor(real_y), preds.cpu(), train_nodes, test_nodes)
print(acc)
print(conf_mat)
print(conf_mat.diagonal()/conf_mat.sum(axis=1))

#pd.Series(preds.cpu().numpy()).to_csv("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/sharp_preds.csv")
conf_labels = pd.Series(sharp.confident_labels)
#conf_labels.to_csv("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/confident_labels.csv")
#int_df = sharp.run_interpretation()
#int_df.to_csv("pbmc_test_interpretation.csv")

