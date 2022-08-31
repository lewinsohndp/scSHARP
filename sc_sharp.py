import utilities
import numpy as np
import pandas as pd
import torch
import os
from gcn_model import GCNModel
import torch

class scSHARP:
    """Class that runs predictions"""

    def __init__(self, data_path, tool_preds, tool_list, marker_path, neighbors=2, config="2_40.txt", ncells="all"):
        self.data_path = data_path
        self.preds_path = tool_preds
        self.tools = tool_list
        self.marker_path = marker_path
        self.neighbors = neighbors
        self.config = config
        self.ncells = ncells
        self.cell_names = None
        self.model = None
        self.final_preds = None
        self.genes = None
        self.X = None
        self.pca_obj = None
        self.batch_size = None
    
    def run_prediction(self, training_epochs=150, thresh=0.51, batch_size=40, seed=8):
        """Trains GCN modle on consensus labels and returns predictions"""
        self.batch_size = batch_size 
        if os.path.exists(self.preds_path):
            all_labels = pd.read_csv(self.preds_path, index_col=0)
            if all_labels.shape[1] != len(self.tools): 
                all_labels = all_labels[self.tools]
                
        else:
            raise Exception("Prediction Dataframe not Found at " + self.data_path)
        
        # read in dataset
        if self.ncells == "all":
            counts = pd.read_csv(self.data_path, index_col=0)
        else:
            counts = pd.read_csv(self.data_path, index_col=0, nrows=self.ncells)
            all_labels = all_labels.head(self.ncells)
        self.X, keep_cells,keep_genes,self.pca_obj = utilities.preprocess(np.array(counts), scale=False, comps=500)
        self.genes = counts.columns.to_numpy()[keep_genes]
        all_labels = all_labels.loc[keep_cells,:]

        _,marker_names = utilities.read_marker_file(self.marker_path)

        self.cell_names = marker_names.copy()
        self.cell_names.sort()

        all_labels_factored = utilities.factorize_df(all_labels, marker_names)
        encoded_labels = utilities.encode_predictions(all_labels_factored)

        confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = thresh)

        train_nodes = np.where(confident_labels != -1)[0]
        test_nodes = np.where(confident_labels == -1)[0]

        dataset  = torch.utils.data.TensorDataset(torch.tensor(self.X), torch.tensor(confident_labels))
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        test_dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)

        if self.model == None: self.model = GCNModel(self.config, neighbors=self.neighbors, target_types=len(marker_names), seed=seed)
        self.model.train(dataloader, epochs=training_epochs, verbose=True)

        preds,_ = self.model.predict(test_dataloader)
        self.final_preds = preds.max(dim=1)[1]

        return self.final_preds, train_nodes, test_nodes, keep_cells

    def run_interpretation(self):
        """Runs model gradient based interpretation"""
        
        int_df = utilities.run_interpretation(self.model, self.X, self.pca_obj, self.final_preds, self.genes, self.batch_size)
        int_df.columns = self.cell_names
        
        return int_df

    def save_model(self, file_path):
        """Save model as serialized object at specified path"""

        torch.save(self.model, file_path)

    def load_model(self, file_path):
        """Load model as serialized object at specified path"""

        self.model = torch.load(file_path)