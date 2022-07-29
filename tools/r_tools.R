# script that runs R prediction tools

require(Seurat)
require(SCINA)
require(scSorter)
require(dplyr)
require(SingleR)

run_scina <- function(data, markers=NULL, ref=NULL){
  if(is.null(markers)){
    return(F)
  }
  
  results <- SCINA(as.matrix(data@assays$RNA@data), markers)
  scina_preds <- results$cell_labels
  
  return(scina_preds)
  
}

run_scsorter <- function(data, markers=NULL, ref=NULL){
  #need top variable genes
  if(is.null(markers)){
    return(F)
  }
  
  types <- list()
  for(i in 1:length(names(markers))){
    types <- append(types, rep(names(markers)[i], length(unlist(markers[i]))))
  }
  
  types <- unlist(types) 
  markers <- unname(unlist(markers))
  anno <- data.frame(Type=types, Marker=markers)
  
  topgenes <- head(VariableFeatures(data), 2000)
  picked_genes <- unique(c(topgenes, anno$Marker))
  expr <- as.matrix(data@assays$RNA@data)
  expr <- expr[rownames(expr) %in% picked_genes,]
  
  rts <- scSorter(expr, anno)
  scsort_preds = rts$Pred_Type
  
  return(scsort_preds)
}

run_sctype <- function(data, markers=NULL, ref=NULL){
  if(is.null(markers)){
    return(F)
  }
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  es.max = sctype_score(scRNAseqData = as.matrix(data@assays$RNA@data), scaled = F, 
                        gs = markers, gs2 = NULL, gene_names_to_uppercase = F) 
  cL_resutls = do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "NA"
  
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    data@meta.data$customclassif[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  sctype_preds = data$customclassif
  
  return(sctype_preds)
}

run_singler <- function(data, ref, ref_labels){
  results <- SingleR(as.matrix(data@assays$RNA@data), as.matrix(ref@assays$RNA@data), ref_labels[,1])
  return(results$pruned.labels)
}

run <- function(data_path, tools, markers=NULL, marker_names=NULL, ref_path=NULL, ref_labels_path=NULL){
  # add ability to take in seurat counts object
  
  counts <- read.csv(data_path, header=T, row.names = 1)
  data <- CreateSeuratObject(t(counts))
  data <- NormalizeData(data)
  data <- ScaleData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- RunPCA(data)
  data <- FindNeighbors(data)
  data <- FindClusters(data)
  
  for(i in 1:length(markers)){
    markers[i] <- list(unlist(markers[i]))
  }
  
  names(markers) <- marker_names
  
  tools <- unlist(tools)
  #markers <- list(Group1=c("Gene1","Gene3"), Group2=c("Gene2"))
  #print(markers)
  results_df <- data.frame(start=rep(0,nrow(counts)))
  
  if("scina" %in% tools){
    results_df$scina <- run_scina(data, markers)
  }
  if("scsorter" %in% tools){
    results_df$scsorter <- run_scsorter(data, markers)
  }
  if("sctype" %in% tools){
    results_df$sctype <- run_sctype(data, markers)
  }
  
  if(!is.null(ref_path) & ("singler" %in% tools)){
    ref_counts <- read.csv(ref_path, header=T, row.names = 1)
    ref_labels = read.csv(ref_labels_path, header=T, row.names=1)
    ref <- CreateSeuratObject(t(ref_counts))
    ref <- NormalizeData(ref)
    results_df$singler <- run_singler(data, ref, ref_labels)
  }
  
  results_df = results_df[-1]
  return(results_df)
}

#data_path <- "~/desktop/conradLab/thesis/scSHARP/simulations/splat_0.4_de/counts.csv"
#results <- run(data_path)

