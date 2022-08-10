import sys
import numpy as np
import pandas as pd
import json
import torch
import utilities
from gcn_model import GCNModel
import os
import statistics

"""script that takes in config file and returns metrics"""

data_folder = sys.argv[1]
output_folder = sys.argv[2]
config_file = "grid_search_files/" + sys.argv[3]
out_file = sys.argv[3]

with open(config_file) as f:
    config = json.load(f)
    f.close()

gcn_file = "configs/" + config['config']
dropout = float(config['dropout'])
batch_size = int(config['batch_size'])
neighbors = int(config['neighbors'])
tools = ["sctype","scsorter","scina","singler", "scpred"]
votes = 3
training_epochs = 150
data_path = data_folder + "query_counts.csv"
ref_path = data_folder + "ref_counts.csv"
ref_label_path = data_folder + "ref_labels.csv"
marker_path = data_folder + "markers.txt"
if os.path.exists(data_folder + "preds.csv"):
    all_labels = pd.read_csv(data_folder + "preds.csv", index_col=0)
    if all_labels.shape[1] != len(tools): raise Exception("wrong amount of tools in file")
else:
    all_labels = utilities.label_counts(data_path,tools,ref_path,ref_label_path,marker_path)

# read in dataset
X = pd.read_csv(data_path, index_col=0)
X, keep_cells = utilities.preprocess(np.array(X), scale=False)

all_labels = all_labels.loc[keep_cells,:]

_,marker_names = utilities.read_marker_file(marker_path)

all_labels_factored = utilities.factorize_df(all_labels, marker_names)
encoded_labels = utilities.encode_predictions(all_labels_factored)

meta_path = data_folder + "query_meta.csv"
metadata = pd.read_csv(meta_path, index_col=0)
real_y = pd.factorize(metadata['Group'], sort=True)[0]
real_y = real_y[keep_cells]

confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = votes)
train_nodes = np.where(confident_labels != -1)[0]
test_nodes = np.where(confident_labels == -1)[0]

dataset  = torch.utils.data.TensorDataset(torch.tensor(X), torch.tensor(confident_labels))
dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)

test_dataset  = torch.utils.data.TensorDataset(torch.tensor(X), torch.tensor(real_y))
test_dataloader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

random_inits = 10
train_accuracies = [0]*random_inits
test_accuracies = [0]*random_inits
total_accuracies = [0]*random_inits

for i in range(random_inits):
    m = GCNModel(gcn_file, neighbors, dropout=dropout)
    m.train(dataloader, training_epochs, verbose=False)
    metrics = m.validation_metrics(test_dataloader, train_nodes, test_nodes)
    total_accuracies[i] = metrics[0]
    train_accuracies[i] = metrics[2]
    test_accuracies[i] = metrics[4]

avg_total = str(statistics.mean(total_accuracies))
avg_train = str(statistics.mean(train_accuracies))
avg_test = str(statistics.mean(test_accuracies))

sd_total = str(statistics.stdev(total_accuracies))
sd_train = str(statistics.stdev(train_accuracies))
sd_test = str(statistics.stdev(test_accuracies))

with open(output_folder + out_file, "w") as output:
    index = out_file.split(".")[0]
    output.write(index + "," + avg_total + "," + avg_train + "," + avg_test + "," + sd_total + "," + sd_train + "," + sd_test + "," + config['config'] + "," + str(config['dropout']) + "," + str(config['batch_size']) +  "," + str(config["neighbors"]) + "\n")
    output.close()

