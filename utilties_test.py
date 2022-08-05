from ctypes import util
import unittest
import pandas as pd
import numpy as np
import utilities
import torch

class TestUtilities(unittest.TestCase):

    def setUp(self):
        self.pred_df = pd.DataFrame({"tool1":["type2","type2", "type1"], "tool2":["type1", "type2", "type3"], "tool3":[pd.NA, "type3", "type3"]})

    def test_marker_parser(self):
        correct_markers = [["Gene13796", "Gene19442", "Gene28124", "Gene22910", "Gene31354"],["Gene19098", "Gene17789", "Gene4696", "Gene25466", "Gene24205"], ["Gene1957", "Gene6552", "Gene931", "Gene25528", "Gene2460"], ["Gene1682", "Gene20645", "Gene17237", "Gene11077", "Gene32186"]]
        correct_marker_names = ["Group1", "Group2", "Group3", "Group4"]

        markers, marker_names = utilities.read_marker_file("simulations/splat_0.5_de_rq/markers_test.txt")

        assert markers == correct_markers
        assert marker_names == correct_marker_names
    
    def test_factorize(self):
        correct_factored = pd.DataFrame({"tool1":[1,1,0], "tool2":[0,1,2], "tool3":[-1, 2, 2]})

        returned = utilities.factorize_df(self.pred_df, ["type1", "type2", "type3"])
        assert correct_factored.equals(returned) == True
    
    def test_encode(self):
        correct_encoded = np.array([[1,1,0],[0,2,1],[1,0,2]])
        returned = utilities.factorize_df(self.pred_df, ["type1", "type2", "type3"])
        encoded = utilities.encode_predictions(returned)

        assert np.array_equal(correct_encoded, encoded) == True
    
    def test_accuracy(self):
        preds = torch.tensor([0,1,1,2])
        real = torch.tensor([0,1,2,2])
        real_accuracy = .75
       
        assert utilities.pred_accuracy(preds, real) == real_accuracy

    def test_consensu(self):
        encoded = np.array([[1,1,0],[0,2,1],[1,0,2]])
        correct_labels = np.array([-1, 1, 2])
        test_labels = utilities.get_consensus_labels(encoded, 2)
        
        assert np.array_equal(test_labels, correct_labels) == True

if __name__ == '__main__':
    unittest.main() 
        
