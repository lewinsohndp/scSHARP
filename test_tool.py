import unittest
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

class TestTools(unittest.TestCase):
    def setUp(self):
        self.markers = [["Gene14207","Gene19034","Gene2396","Gene948"],["Gene9155","Gene5388","Gene4361","Gene10982"],["Gene13075","Gene142","Gene3523"],["Gene16021","Gene9929","Gene11519","Gene6797"]]
        self.marker_names = ["Group1", "Group2","Group3", "Group4"]
        self.query_path = "simulations/splat_0.2_de_rq/query_counts.csv"
        self.ref_path = "simulations/splat_0.2_de_rq/ref_counts.csv"
        self.ref_label_path = "simulations/splat_0.2_de_rq/ref_labels.csv"

    def test_rtools(self):
        ro.r.source('tools/r_tools.R')
        preds = ro.r.run(self.query_path, ["scsorter", "sctype", "scina", "singler", "scpred"], self.markers, self.marker_names, self.ref_path, self.ref_label_path)
        with localconverter(ro.default_converter + pandas2ri.converter):
            preds = ro.conversion.rpy2py(preds)

        print(preds.head)
        assert preds.shape == (1000,5)
    
    def test_pytools(self):
        pass

if __name__ == '__main__':
    unittest.main() 
        
