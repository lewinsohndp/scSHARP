library(Seurat)

# read in data path as arg 1
args <- commandArgs(trailingOnly = T)
data_path <- args[1]
output_path <- args[2]

counts <- read.csv(data_path, header=T, row.names = 1)
data <- CreateSeuratObject(t(counts), min.cells=3, min.features=200)
data <- NormalizeData(data)
data <- ScaleData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- RunPCA(data)
data <- FindNeighbors(data)
data <- FindClusters(data)
data <- RunUMAP(data, dims = 1:50)

saveRDS(data, file = paste(output_path, "process_seur.rds", sep=""))
write.csv(Embeddings(data, reduction="umap"), paste(output_path,"umap_embedding.csv",sep=""))
