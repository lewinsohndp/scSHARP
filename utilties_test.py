from ctypes import util
import unittest
import pandas as pd
import utilities

class TestUtilities(unittest.TestCase):

    def test_marker_parser(self):
        correct_markers = [["Gene14207", "Gene19034", "Gene2396", "Gene948", "Gene9141"],["Gene9155", "Gene5388", "Gene4361", "Gene10982", "Gene16515"], ["Gene13075", "Gene142", "Gene3523", "Gene14279", "Gene14643"], ["Gene16021", "Gene9929", "Gene11519", "Gene6797", "Gene13413"]]
        correct_marker_names = ["Group1", "Group2", "Group3", "Group4"]

        markers, marker_names = utilities.read_marker_file("simulations/splat_0.2_de_rq/markers.txt")

        assert markers == correct_markers
        assert marker_names == correct_marker_names

if __name__ == '__main__':
    unittest.main() 
        
