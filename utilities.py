import torch
from torch_geometric.nn import MessagePassing
from torch.nn import Sequential as Seq, Linear, SiLU, Dropout, ELU, Tanh
import anndata as ad
import scanpy as sc
from sklearn.decomposition import PCA
import numpy as np
import random
import math
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import json
import pandas as pd
from sklearn.preprocessing import OneHotEncoder
from captum.attr import IntegratedGradients, DeepLift
import math

"""General functions and definitions"""

class EdgeConv(MessagePassing):
    """Edge convolutional layer definition"""

    def __init__(self, in_channels, out_channels):
        super().__init__(aggr='add') #  "Max" aggregation.
        self.mlp = Seq(Linear(2 * in_channels, out_channels),
                       SiLU(),
                       Linear(out_channels, out_channels))

    def forward(self, x, edge_index):
        # x has shape [N, in_channels]
        # edge_index has shape [2, E]

        return self.propagate(edge_index, x=x)

    def message(self, x_i, x_j):
        # x_i has shape [E, in_channels]
        # x_j has shape [E, in_channels]

        tmp = torch.cat([x_i, x_j - x_i], dim=1)  # tmp has shape [E, 2 * in_channels]
        return self.mlp(tmp)

def preprocess(data, normalize=True, scale=True, targetsum=1e4, run_pca=True, comps=500):
    """method for preprocessing raw counts matrix"""

    adata = ad.AnnData(data, dtype=data.dtype)
    row_filter, _ = sc.pp.filter_cells(adata, min_genes=200, inplace=False)
    col_filter,_ = sc.pp.filter_genes(adata, min_cells=3, inplace=False)
    subset = adata[row_filter,col_filter]
    adata = ad.AnnData(subset.X, dtype=subset.X.dtype)
    if normalize:
        sc.pp.normalize_total(adata, target_sum=targetsum, layer=None)
    new_data = adata.X
    if normalize:
        new_data = sc.pp.log1p(new_data)
    
    if scale:
        new_data = sc.pp.scale(new_data)
    
    pca = None
    if run_pca:
        pca = PCA(n_components=comps, random_state=8)
        new_data = pca.fit_transform(new_data)

    return new_data, row_filter, col_filter, pca

def mask_labels(labels, masking_pct):
    """method for masking labels"""
    random.seed(8)
    subset = random.sample(range(len(labels)), math.floor((len(labels)*masking_pct)))
    masked_labels = np.copy(labels)
    masked_labels[subset] = -1 

    return labels, masked_labels  

def read_marker_file(file_path):
    """parse marker file"""

    marker_df = pd.read_csv(file_path, header=None, index_col=0)
    markers = np.array(marker_df).tolist()
    marker_names = list(marker_df.index)
    """with open(file_path) as f:
        lines = f.readlines()
        f.close()
    
    markers = [[] for _ in range(len(lines))]
    marker_names = [""]*len(markers)
    for i, line in enumerate(lines):
        temp_name = line.split(sep=":")[0]
        temp_markers = line.split(sep=":")[1].split(sep=",")
        temp_markers = [s.strip() for s in temp_markers]
        markers[i] = temp_markers
        marker_names[i] = temp_name"""

    return markers, marker_names

def label_counts(data_path, tools, ref_path, ref_label_path, marker_path):
    """Label inputted dataset with mutlitple annotation tools"""

    markers, marker_names = read_marker_file(marker_path)

    ro.r.source('tools/r_tools.R')
    #markers = 
    preds = ro.r.run(data_path, tools, markers, marker_names, ref_path, ref_label_path)
    with localconverter(ro.default_converter + pandas2ri.converter):
        preds = ro.conversion.rpy2py(preds)
    
    return preds

def load_model(file_path):
        """loads model from json format"""
        
        f = open(file_path)
        data = json.load(f)
        f.close()
        final_layers = []

        for layer in data['layers']:
            if layer['type'] == 'EdgeConv':
                new_layer = EdgeConv(layer['input'],layer['output'])
            
            elif layer['type'] == 'Linear':
                new_layer = torch.nn.Linear(layer['input'], layer['output'])
            
            elif layer['type'] == 'Sigmoid':
                new_layer = torch.nn.Sigmoid()
            
            elif layer['type'] == 'BatchNorm':
                new_layer = torch.nn.BatchNorm1d(layer['input'])

            else:
                raise Exception("Unrecognizable layer type:" + layer['type'])
            
            final_layers.append(new_layer)
        
        return final_layers

def factorize_df(df, all_cells):
    """factorizes all columns in pandas df"""

    # add array with all cell options so factor is the same for each column
    dummy = []
    for i in range(len(all_cells)):
        dummy.append([all_cells[i]]*df.shape[1])
    dummy = pd.DataFrame(dummy)
    dummy.columns = df.columns
    df = pd.concat([df,dummy])
    
    temp = df.apply(pd.factorize, axis=0, sort=True)
    temp = temp.iloc[0,:]
    indices = list(temp.index)
    d = {key: None for key in indices}
    for i in range(temp.shape[0]):
        d[indices[i]] = temp.iloc[i]

    return pd.DataFrame(d).iloc[:-len(all_cells)]

def encode_predictions(df):
    """encodes predictions for each cell with 1 for each prediction"""
    all_preds = []
    for i in range(df.shape[1]):
        all_preds.append(df.iloc[:,i].to_numpy())
    all_preds = np.array(all_preds).flatten()
    
    #add -1 then remove so encoder takes into account unknowns even if there isn't any
    all_preds = np.append(all_preds, -1)
    enc = OneHotEncoder(drop='first')
    encoded_y = enc.fit_transform(all_preds.reshape(-1,1)).toarray()
    encoded_y = encoded_y[:-1,:]
    # need to add three scores together
    final_encoded = np.zeros(shape=(df.shape[0],encoded_y.shape[1]))
    scoring_length = df.shape[0]
    lower =0
    upper = scoring_length
    for i in range(int(len(encoded_y)/df.shape[0])):
        final_encoded += encoded_y[lower:upper,:]
        lower = upper
        upper += scoring_length
    
    return final_encoded

def pred_accuracy(preds, real):
    """returns accuracy of predictions"""
    return float((torch.tensor(preds) == torch.tensor(real)).type(torch.FloatTensor).mean().numpy())

def get_consensus_labels(encoded_y, necessary_vote):
    """method that gets consensus vote of multiple prediction tools
    If vote is < 1 then taken as threshold pct to be >= to
    """
    confident_labels = np.zeros(shape = (encoded_y.shape[0],))
    for i in range(encoded_y.shape[0]):
        row = encoded_y[i,:]

        if necessary_vote < 1:
            row = row / row.sum()
        
        max_index = np.argmax(row)
        if row[max_index] >= (necessary_vote):
            confident_labels[i] = max_index
        else: confident_labels[i] = -1
    
    return confident_labels

def filter_scores(scores, thresh = 0.5):
    """filters out score columns with NAs > threshold"""
    keep_cols = []
    for col in scores.columns:
        pct_na = (scores[col].isna().sum()) / scores.shape[0]

        if pct_na < thresh: keep_cols.append(col)
    
    return scores[keep_cols]

def run_interpretation(model, X, pca_obj, predictions, genes):
    """Method to run interpretation on model"""

    predictions = predictions.cpu()
    prediction_names = predictions.unique().tolist()
    classes = [None]*len(prediction_names)
    for i, pred in enumerate(prediction_names):
        classes[i] = np.where(predictions == pred)[0]
    
    # fix adding eval mode

    dl = DeepLift(model)
    baseline = torch.FloatTensor(np.full(X.shape, X.min()))

    attributions = np.zeros((len(prediction_names), X.shape[0], X.shape[1]))

    for i, pred_name in enumerate(prediction_names):
        attributions[i] = dl.attribute(torch.FloatTensor(X).to(model.device), baseline.to(model.device), target=pred_name, return_convergence_delta=True)[0].cpu().detach()

    mean_attributions = np.zeros((len(prediction_names), X.shape[1]))
    for i in range(attributions.shape[0]):
        mean_attributions[i] = torch.mean(torch.tensor(attributions[i]), 0)
    
    mean_attributions = torch.FloatTensor(mean_attributions.T)

    cor_loadings = torch.FloatTensor(pca_obj.components_.T * np.sqrt(pca_obj.explained_variance_))
    gene_att_scores = torch.mm(cor_loadings, mean_attributions)
    gene_att_scores = gene_att_scores.numpy()
    att_df = pd.DataFrame(gene_att_scores)
    att_df.index = genes
    att_df.columns = prediction_names

    return att_df
