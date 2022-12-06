from . import utilities
import torch

class PCAModel(torch.nn.Module):
    """class for pca model that takes in pre fitted components"""

    def __init__(self, pca_comps, pca_mean):
        super().__init__()
        use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda:0" if use_cuda else "cpu")
        self.pca_comps = torch.nn.Parameter(torch.FloatTensor(pca_comps).to(self.device))
        self.pca_mean = torch.nn.Parameter(torch.FloatTensor(pca_mean).to(self.device))
        use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda:0" if use_cuda else "cpu")

    def forward(self, X):
        """forward pass for pca transformation"""
        #X = torch.FloatTensor(X).to(self.device)
        return torch.mm(X-self.pca_mean, self.pca_comps.T)
