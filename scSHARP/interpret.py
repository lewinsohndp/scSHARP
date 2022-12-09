import torch
import pandas as pd
import numpy as np
from captum.attr import IntegratedGradients, DeepLift, DeepLiftShap, FeaturePermutation

# TODO check with Daniel that we can kill the old interpret model method or at least rename it

def interpret_model_new(model, X, predictions, genes, batch_size, device, batches=None):
    dataset  = torch.utils.data.TensorDataset(torch.FloatTensor(X), predictions.cpu())
    dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    predictions = predictions.cpu()
    prediction_names = predictions.unique().tolist()
    classes = [None]*len(prediction_names)
    for i, pred in enumerate(prediction_names):
        classes[i] = np.where(predictions == pred)[0]

    dl = DeepLift(model)
    #dl = FeaturePermutation(model)
    #attributions = np.zeros((X.shape[0], X.shape[1]))
    temp_atts = None
    counter = 0
    for data,preds in dataloader:
        baseline = torch.FloatTensor(np.zeros(data.shape))
        temp = dl.attribute(data.to(device), baseline.to(device), target=preds.to(device), return_convergence_delta=True)[0].cpu().detach()
        #temp = dl.attribute(data.to(device), target=preds.to(device)).cpu().detach()
           
        if temp_atts == None: temp_atts = temp
        else:
            temp_atts = torch.cat((temp_atts, temp), 0)
        #if counter % 10 == 0: print(str(counter))
        counter +=1 
        if batches != None:
            if counter == batches: break
    attributions = temp_atts
    
    mean_attributions = np.zeros((len(prediction_names), X.shape[1]))
    for i in range(len(prediction_names)):
        mean_attributions[i] = torch.mean(attributions[classes[i][classes[i] < attributions.shape[0]],:], 0)
    
    mean_attributions = mean_attributions.T
    att_df = pd.DataFrame(mean_attributions)
    att_df.index = genes
    att_df.columns = prediction_names

    return att_df

def heat_map(att_df):
    
    pass

def interpret_model(model, X, pca_obj, predictions, genes, batch_size):
    """Method to run interpretation on model"""
    
    dataset  = torch.utils.data.TensorDataset(torch.FloatTensor(X), predictions.cpu())
    dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    predictions = predictions.cpu()
    prediction_names = predictions.unique().tolist()
    classes = [None]*len(prediction_names)
    for i, pred in enumerate(prediction_names):
        classes[i] = np.where(predictions == pred)[0]
    #pd.DataFrame(classes).to_csv("preds_test.csv")
    # fix adding eval mode

    #dl = DeepLift(model)
    #dl = DeepLiftShap(model)
    dl = FeaturePermutation(model)
    baseline = torch.FloatTensor(np.full(X.shape, X.min()))

    #attributions = np.zeros((len(prediction_names), X.shape[0], X.shape[1]))
    attributions = np.zeros((X.shape[0], X.shape[1]))
    temp_atts = None
    for data,preds in dataloader:
        baseline = torch.FloatTensor(np.full(data.shape, X.min()))
        #baseline = torch.FloatTensor(np.zeros(data.shape))
        #baseline = torch.FloatTensor(np.full(data.shape, X.max()))
        #temp = dl.attribute(data.to(model.device), baseline.to(model.device), target=preds.to(model.device), return_convergence_delta=True)[0].cpu().detach()
        temp = dl.attribute(data.to(model.device), target=preds.to(model.device)).cpu().detach()
        if temp_atts == None: temp_atts = temp
        else:
            temp_atts = torch.cat((temp_atts, temp), 0)
    
    attributions = temp_atts
    #pd.DataFrame(attributions.numpy()).to_csv("single_atts.csv")
    mean_attributions = np.zeros((len(prediction_names), X.shape[1]))
    for i in range(len(prediction_names)):
        mean_attributions[i] = torch.mean(attributions[classes[i],:], 0)
    
    mean_attributions = torch.FloatTensor(mean_attributions.T)

    cor_loadings = torch.FloatTensor(pca_obj.components_.T * np.sqrt(pca_obj.explained_variance_))
    #gene_att_scores = torch.mm(cor_loadings, mean_attributions)
    gene_att_scores = torch.mm(torch.FloatTensor(pca_obj.components_.T), mean_attributions)
    gene_att_scores = gene_att_scores.numpy()
    att_df = pd.DataFrame(gene_att_scores)
    att_df.index = genes
    att_df.columns = prediction_names
    
     
    #pd.DataFrame(cor_loadings.numpy()).to_csv("cor_loadings.csv")
    #pd.DataFrame(mean_attributions.numpy()).to_csv("mean_atts.csv")
    #pd.DataFrame(pca_obj.components_.T).to_csv("comp_loadings.csv")
    return att_df