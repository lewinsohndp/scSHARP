library(splatter)
library(Seurat)
library(scater)
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

folder <- "~/desktop/conradLab/thesis/scSHARP/simulations/"
for(i in seq(.1,.2,length.out=2)){
  full_path <- paste(folder,"splat_",i,"_de_rq/", sep="")
  dir.create(full_path)
  
  sim.groups <- splatSimulate(batchCells=c(1000, 1000), group.prob = c(.423, .227, .226, .124), method = "groups", nGenes = 20000, de.facScale=i, verbose = FALSE, seed=8)
  
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
  umap <- DimPlot(whole, reduction="umap")
  ggsave(paste(full_path,"query_umap.pdf",sep=""), plot=umap, device="pdf")
  umap <- DimPlot(whole, reduction="umap", group.by = "Group")
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
  
  write.csv(t(whole@assays$RNA@counts), paste(full_path,"whole_counts.csv",sep=""))
  write.csv(whole@meta.data, paste(full_path,"whole_meta.csv",sep=""))
  
  write.csv(t(ref@assays$RNA@counts), paste(full_path,"ref_counts.csv",sep=""))
  write.csv(ref@meta.data, paste(full_path,"ref_meta.csv",sep=""))
  write.csv(data.frame(ref$Group), paste(full_path, "ref_labels.csv", sep=""))
  
  write.csv(t(query@assays$RNA@counts), paste(full_path,"query_counts.csv",sep=""))
  write.csv(query@meta.data, paste(full_path,"query_meta.csv",sep=""))
}
