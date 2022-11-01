library(splatter)
library(Seurat)
#library(scater)
library(ggplot2)

run_analysis <- function(seur_obj){
  seur_obj <- NormalizeData(seur_obj)
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)
  seur_obj <- ScaleData(seur_obj)
  seur_obj <- RunPCA(seur_obj)
  seur_obj <- FindNeighbors(seur_obj)
  seur_obj <- FindClusters(seur_obj)
  seur_obj <- RunUMAP(seur_obj, dims = 1:50)
  return(seur_obj)
}

folder <- "/home/groups/ConradLab/daniel/sharp_data/sharp_sims_de_testing/"
param_data_path <- "~/scSHARP/filtered_gene_bc_matrices/GRCh38/matrix.mtx"
# folder <- "~/desktop/conradlab/thesis/scsharp/simulations/"
# param_data_path <- "~/desktop/conradlab/thesis/scSHARP/filtered_gene_bc_matrices/GRCh38/matrix.mtx"

counts <- Matrix::readMM(param_data_path)
set.seed(8)
# counts <- counts[,sample(1:ncol(counts), 1000)]
# filtering for acutally expressed genes
# counts = counts[rowSums(counts)>1, ]

params <- splatEstimate(as.matrix(counts))

for(i in seq(.3,1.0,length.out=8)){
  full_path <- paste(folder,"splat_",i,"_de_rq_103122/", sep="")
  dir.create(full_path)
  params <- setParams(params, dropout.type='experiment', batchCells=c(1000,1000), batch.facScale=0.5, group.prob = c(.25,.25,.25,.25), de.facScale=i, seed=888)
  print(params)
  sim.groups <- splatSimulate(params, method = "groups", verbose = T)
  
  whole <- CreateSeuratObject(assays(sim.groups)$counts, meta.data=as.data.frame(sim.groups@colData))
  ref <- subset(x = whole, subset = Batch == "Batch1")
  query <- subset(x = whole, subset = Batch == "Batch2")
  
  whole <- run_analysis(whole)
  ref <- run_analysis(ref)
  query <- run_analysis(query)
  
  #whole umaps
  umap <- DimPlot(whole, reduction="umap")
  ggsave(paste(full_path,"whole_umap.pdf",sep=""), plot=umap, device="pdf")
  umap <- DimPlot(whole, reduction="umap", group.by = "Group")
  ggsave(paste(full_path,"whole_umap_groups.pdf",sep=""), plot=umap, device="pdf")
  umap <- DimPlot(whole, reduction="umap", group.by = "Batch")
  ggsave(paste(full_path,"whole_umap_batches.pdf",sep=""), plot=umap, device="pdf")
  
  #ref umaps
  umap <- DimPlot(ref, reduction="umap")
  ggsave(paste(full_path,"ref_umap.pdf",sep=""), plot=umap, device="pdf")
  umap <- DimPlot(ref, reduction="umap", group.by = "Group")
  ggsave(paste(full_path,"ref_umap_groups.pdf",sep=""), plot=umap, device="pdf")
  
  #query umaps
  umap <- DimPlot(query, reduction="umap")
  ggsave(paste(full_path,"query_umap.pdf",sep=""), plot=umap, device="pdf")
  umap <- DimPlot(query, reduction="umap", group.by = "Group")
  ggsave(paste(full_path,"query_umap_groups.pdf",sep=""), plot=umap, device="pdf")
  
  Idents(object = whole) <- "Group"
  group1.markers <- FindMarkers(whole, ident.1 = "Group1", ident.2 = NULL, only.pos = TRUE)[1:20,]
  group2.markers <- FindMarkers(whole, ident.1 = "Group2", ident.2 = NULL, only.pos = TRUE)[1:20,]
  group3.markers <- FindMarkers(whole, ident.1 = "Group3", ident.2 = NULL, only.pos = TRUE)[1:20,]
  group4.markers <- FindMarkers(whole, ident.1 = "Group4", ident.2 = NULL, only.pos = TRUE)[1:20,]
  write.csv(group1.markers, paste(full_path,"group1_de.csv",sep=""))
  write.csv(group2.markers, paste(full_path,"group2_de.csv",sep=""))
  write.csv(group3.markers, paste(full_path,"group3_de.csv",sep=""))
  write.csv(group4.markers, paste(full_path,"group4_de.csv",sep=""))
  
  # write markers.txt
  set.seed(8)
  group1 <- paste("Group1,",paste(row.names(group1.markers)[sample(1:10, 5)], collapse=","), sep="")
  group2 <- paste("Group2,",paste(row.names(group2.markers)[sample(1:10, 5)], collapse=","), sep="")
  group3 <- paste("Group3,",paste(row.names(group3.markers)[sample(1:10, 5)], collapse=","), sep="")
  group4 <- paste("Group4,",paste(row.names(group4.markers)[sample(1:10, 5)], collapse=","), sep="")
  all <- c(group1, group2, group3, group4)
  write(all, paste(full_path, "markers.txt", sep=""), ncolumns = 1, sep = "\n")
  
  write.csv(t(whole@assays$RNA@counts), paste(full_path,"whole_counts.csv",sep=""))
  write.csv(whole@meta.data, paste(full_path,"whole_meta.csv",sep=""))
  
  write.csv(t(ref@assays$RNA@counts), paste(full_path,"ref_counts.csv",sep=""))
  write.csv(ref@meta.data, paste(full_path,"ref_meta.csv",sep=""))
  write.csv(data.frame(ref$Group), paste(full_path, "ref_labels.csv", sep=""))
  
  write.csv(t(query@assays$RNA@counts), paste(full_path,"query_counts.csv",sep=""))
  write.csv(query@meta.data, paste(full_path,"query_meta.csv",sep=""))
}
