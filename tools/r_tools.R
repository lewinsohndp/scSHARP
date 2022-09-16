# script that runs R prediction tools

run_scina <- function(data, markers=NULL, ref=NULL){
  # can return "unknown"
  
  if(is.null(markers)){
    return(F)
  }
  print("Running SCINA") 
  norm_counts <- GetAssayData(object = data, slot = "data") 
  results <- SCINA(norm_counts, markers, rm_overlap=F)
  scina_preds <- results$cell_labels
  
  scina_preds = replace(scina_preds, scina_preds=="unknown", NA)
  return(scina_preds)
  
}

run_scsorter <- function(data, markers=NULL, ref=NULL){
  #need top variable genes
  #can return "Unknown"
  
  if(is.null(markers)){
    return(F)
  }
  
  print("Running scSorter") 
  types <- list()
  for(i in 1:length(names(markers))){
    types <- append(types, rep(names(markers)[i], length(unlist(markers[i]))))
  }
  
  types <- unlist(types) 
  markers <- unname(unlist(markers))
  anno <- data.frame(Type=types, Marker=markers)
  
  topgenes <- head(VariableFeatures(data), 2000)
  picked_genes <- unique(c(topgenes, anno$Marker))
  #expr <- as.matrix(data@assays$RNA@data)
  expr <- as.matrix(GetAssayData(object = data, slot = "data"))
  expr <- expr[rownames(expr) %in% picked_genes,]
  
  rts <- scSorter(expr, anno)
  scsort_preds = rts$Pred_Type
  scsort_preds = replace(scsort_preds, scsort_preds=="Unknown", NA)
  
  return(scsort_preds)
}

run_sctype <- function(data, markers=NULL, ref=NULL){
  if(is.null(markers)){
    return(F)
  }
  print("Running scType")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  norm_counts <- as.matrix(GetAssayData(object = data, slot = "data"))
  es.max = sctype_score(scRNAseqData = norm_counts, scaled = F, 
                        gs = markers, gs2 = NULL, gene_names_to_uppercase = F) 
  cL_resutls = do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = NA
  
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    data@meta.data$customclassif[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  sctype_preds = data$customclassif
  
  return(sctype_preds)
}

run_singler <- function(data, ref, ref_labels){
  norm_counts <- as.matrix(GetAssayData(object = data, slot = "data"))
  results <- SingleR(norm_counts, as.matrix(ref@assays$RNA@data), ref_labels[,1])
  return(results$pruned.labels)
}

run_scpred <- function(data, ref, ref_labels){
  # can return "unassigned"
  ref <- FindVariableFeatures(ref)
  ref <- ScaleData(ref)
  ref <- RunPCA(ref)
  ref <- RunUMAP(ref, dims=1:30)
  ref <- AddMetaData(ref, ref_labels, col.name="celltype")
  ref <- getFeatureSpace(ref, "celltype")
  ref <- trainModel(ref)
  data <- scPredict(data, ref)
  
  scpreds = data$scpred_prediction
  scpreds = replace(scpreds, scpreds=="unassigned", NA)
  
  return(scpreds)
}

run <- function(data_path, tools, markers=NULL, marker_names=NULL, ref_path=NULL, ref_labels_path=NULL){
  # add ability to take in seurat counts object
  
  require(Seurat)
  require(SCINA)
  require(scSorter)
  require(dplyr)
  require(SingleR)
  require("scPred")
  #set.seed(25) 
  counts <- read.csv(data_path, header=T, row.names = 1)
  print(dim(counts))
  data <- CreateSeuratObject(t(counts), min.cells=3, min.features=200)
  #data <- CreateSeuratObject(t(counts))
  print(dim(data@assays$RNA@counts))
  data <- NormalizeData(data)
  data <- ScaleData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- RunPCA(data)
  data <- FindNeighbors(data)
  data <- FindClusters(data)
  
  markers <- unname(markers)
  for(i in 1:length(markers)){
    markers[i] <- list(unlist(markers[i]))
  }
  
  names(markers) <- marker_names
  
  tools <- unlist(tools)
  #markers <- list(Group1=c("Gene1","Gene3"), Group2=c("Gene2"))
  #print(markers)
  results_df <- data.frame(start=rep(0,ncol(x = data)))
  #print(row.names(data@assays$RNA@data))
  if("scina" %in% tools){
    results_df$scina <- run_scina(data, markers)
  }
  if("scsorter" %in% tools){
    results_df$scsorter <- run_scsorter(data, markers)
  }
  if("sctype" %in% tools){
    results_df$sctype <- run_sctype(data, markers)
  }
  
  if((!is.null(ref_path)) & (!is.na(ref_path))){
    ref_counts <- read.csv(ref_path, header=T, row.names = 1)
    ref_labels = read.csv(ref_labels_path, header=T, row.names=1)
    ref <- CreateSeuratObject(t(ref_counts))
    ref <- NormalizeData(ref)
    
    if("singler" %in% tools){
      results_df$singler <- run_singler(data, ref, ref_labels)
    }
    if("scpred" %in% tools){
      results_df$scpred <- run_scpred(data, ref, ref_labels)
    }
    
  }
  
  results_df = results_df[-1]
  row.names(results_df) <- row.names(data@meta.data)
  return(results_df)
}

# args: data path, out path, tools (comma separated), marker_path, ref_path, ref_label_path
args <- commandArgs(trailingOnly = T)
if(length(args) > 2){
  data_path <- args[1]
  out_path <- args[2]
  tools <- args[3]
  marker_path <- args[4]
  ref_path <- args[5]
  ref_label_path <- args[6]
  
  tools <- unlist(strsplit(tools,","))
  
  print(tools)
  marker_df <- read.csv(marker_path, header=F, row.names = 1)
  markers <- as.list(as.data.frame(t(marker_df)))
  markers <- lapply(markers, function(z){ z[!is.na(z) & z != ""]})
  marker_names <- row.names(marker_df)
  print(markers)
  print(marker_names)
  print(ref_path)
  results <- run(data_path, tools, markers, marker_names, ref_path, ref_label_path)
  
  write.csv(results, paste(out_path,"preds.csv", sep=""))
  
}




