library(Seurat)

# read in data path as arg 1
args <- commandArgs(trailingOnly = T)
data_path <- args[1]

data <- readRDS(paste(data_path,"process_seur.rds",sep=""))
Idents(data) <- "sharp_preds"
print(dim(data))
data.markers <- FindAllMarkers(data, only.pos = TRUE, verbose=TRUE, random.seed=1)
write.csv(data.markers, paste(data_path, "de_genes.csv",sep=""))
