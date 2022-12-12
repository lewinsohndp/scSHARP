from scSHARP.sc_sharp import scSHARP

data_path = "simulations/splat_0.7/query_counts.csv.gz"
tool_preds = "simulations/splat_0.7/preds.csv"
tool_list = ["scina", "scsorter", "sctype", "scpred", "singler"]
out_path = "./output"
marker_path = "simulations/splat_0.7/markers.txt"
neighbors=2
config="configs/2_25.txt"
ref_path = "examples/splat_0.8/ref_counts.csv.gz"
ref_label_path = "examples/splat_0.8/ref_labels.csv"

sharp = scSHARP(data_path, tool_list, marker_path, neighbors, config)

sharp.run_tools(out_path, ref_path, ref_label_path)