# scSHARP

Repository for "Consensus Label Propagation with Graph Convolutional Networks for Single-Cell RNA Sequencing Cell Type Annotation" extended abstract submission for Learning on Graphs Conference and full paper available at https://doi.org/10.1101/2022.11.23.517739.

## To use:
First clone repository.

Next, run tools/r_tools.R to generate predictions from scType, SCINA, scSorter, scPred, and SingleR. Example command:  
Rscript tools/r_tools.R path/to/data_csv path/to/output/ tool1,tool2,tool3 path/to/markers path/to/reference_csv path/to/reference_labels

Reference data set and labels are not required, but if you are not using them, be sure to not specify any reference-based tools in tool list.  
Marker list should be in CSV format:  
cell_type_1,gene1,gene2  
cell_type_2,gene3,gene4  

Following generation of tool predictions we are ready to train the model and create predictions.  
Prepare model with the following code:  
data_path = "path/to/data_csv"  
tool_preds = "path/to/tool/predictions"  
tool_list = ["tool1", "tool2", "tool3"] # tools to use for model. 
marker_path = "path/to/markers"  
neighbors=2 # neighbors to use in graph generation. 
config="configs/2_25.txt" # hyperparameter configuration file to use  
sharp = scSHARP(data_path, tool_preds, tool_list, marker_path, neighbors, config)  

preds, train_nodes, test_nodes, keep_cells = sharp.run_prediction(training_epochs=150, thresh=0.51, batch_size=20, seed=8)  

preds is the list of cell type predictions, train_nodes are the cell designated as confidently labeled, test_nodes are cells designated as unconfidently labeled, and keep_cells is the list of cells not filtered by preprocessing.

The following code runs DeepLIFT interpretation:  
int_df = sharp.run_interpretation()  
int_df = int_df.abs()  
scale_int_df = pd.DataFrame(preprocessing.scale(int_df, with_mean=False))  
scale_int_df.columns = int_df.columns  
scale_int_df.index = int_df.index  

Please note, this pipeline currently does some initial filtering removing genes expressed in fewer than 3 cells and cells expressed in fewer than 200 genes.
