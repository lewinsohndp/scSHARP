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
        correct_markers = [["Gene32664","Gene10522",'Gene15390',"Gene7622","Gene13270"],["Gene15787","Gene7622","Gene14624","Gene3116","Gene21573"],["Gene7733","Gene29807","Gene19010","Gene26085","Gene3915"],["Gene7884","Gene30004","Gene22657","Gene12393","Gene22733"]]
        correct_marker_names = ["Group1", "Group2", "Group3", "Group4"]

        markers, marker_names = utilities.read_marker_file("simulations/splat_0.6_de_rq/markers.txt")

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

    def test_consensus(self):
        encoded = np.array([[1,1,0],[0,2,1],[1,0,2]])
        correct_labels = np.array([-1, 1, 2])
        test_labels = utilities.get_consensus_labels(encoded, 2)
        
        assert np.array_equal(test_labels, correct_labels) == True
    
    def test_consensus_pct(self):
        encoded = np.array([[1,1,0],[0,2,1],[1,0,2]])
        correct_labels = np.array([-1, 1, 2])
        test_labels = utilities.get_consensus_labels(encoded, .51)
        
        assert np.array_equal(test_labels, correct_labels) == True

if __name__ == '__main__':
    unittest.main() 
        
