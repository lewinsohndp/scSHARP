from . import utilities
#import scSHARP.utilities as utilities
import numpy as np
import pandas as pd
import torch
import os
from .gcn_model import GCNModel
import torch
from .pca_model import PCAModel
from sklearn.decomposition import PCA

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
        self.counts = None
        self.keep_cells = None
        self.confident_labels = None
    
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
            self.counts = pd.read_csv(self.data_path, index_col=0)
        else:
            self.counts = pd.read_csv(self.data_path, index_col=0, nrows=self.ncells)
            all_labels = all_labels.head(self.ncells)
        self.X, self.keep_cells,keep_genes,self.pca_obj = utilities.preprocess(np.array(self.counts), scale=False, comps=500)
        self.genes = self.counts.columns.to_numpy()[keep_genes]
        #all_labels = all_labels.loc[self.keep_cells,:]

        _,marker_names = utilities.read_marker_file(self.marker_path)

        self.cell_names = marker_names.copy()
        self.cell_names.sort()

        all_labels_factored = utilities.factorize_df(all_labels, marker_names)
        encoded_labels = utilities.encode_predictions(all_labels_factored)

        self.confident_labels = utilities.get_consensus_labels(encoded_labels, necessary_vote = thresh)

        train_nodes = np.where(self.confident_labels != -1)[0]
        test_nodes = np.where(self.confident_labels == -1)[0]

        dataset  = torch.utils.data.TensorDataset(torch.tensor(self.X), torch.tensor(self.confident_labels))
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        test_dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)

        if self.model == None: self.model = GCNModel(self.config, neighbors=self.neighbors, target_types=len(marker_names), seed=seed)
        self.model.train(dataloader, epochs=training_epochs, verbose=True)

        preds,_ = self.model.predict(test_dataloader)
        self.conf_scores, self.final_preds = preds.max(dim=1)

        return self.final_preds, train_nodes, test_nodes, self.keep_cells, self.conf_scores

    def run_interpretation(self):
        """Runs model gradient based interpretation"""
        X,_,_,_ = utilities.preprocess(np.array(self.counts), scale=False, run_pca=False)
        pca = PCA(n_components=500, random_state=8)
        pca.fit(X)
        pca_mod = PCAModel(pca.components_, pca.mean_)
        seq = torch.nn.Sequential(pca_mod, self.model)
        #meta_path = "/home/groups/ConradLab/daniel/sharp_data/pbmc_test/labels_cd4-8.csv"
        #metadata = pd.read_csv(meta_path, index_col=0)
        #real_y = pd.factorize(metadata.iloc[:,0], sort=True)[0]
        #real_y = real_y[self.keep_cells]
        #real_y = torch.tensor(real_y)
        int_df = utilities.run_interpretation_new(seq, X, self.final_preds, self.genes, self.batch_size, self.model.device)
        int_df.columns = self.cell_names
        
        return int_df

    def save_model(self, file_path):
        """Save model as serialized object at specified path"""

        torch.save(self.model, file_path)

    def load_model(self, file_path):
        """Load model as serialized object at specified path"""

        self.model = torch.load(file_path)
