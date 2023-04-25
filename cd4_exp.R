library(Seurat)

seur_obj <- readRDS("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/process_seur.rds")
real_labels <- read.csv("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/labels_cd4-8.csv", row.names=1)
seur_obj$sharp_preds[seur_obj$sharp_preds == 0] = "b_cells"
seur_obj$sharp_preds[seur_obj$sharp_preds == 1] = "cd14_monocytes"
seur_obj$sharp_preds[seur_obj$sharp_preds == 2] = "cd4_t_cell"
seur_obj$sharp_preds[seur_obj$sharp_preds == 3] = "cd56_nk"
seur_obj$sharp_preds[seur_obj$sharp_preds == 4] = "cd8_t_cell"
seur_obj$sharp_preds <- factor(seur_obj$sharp_preds, levels = c("b_cells", "cd14_monocytes", "cd56_nk", "cd8_t_cell", "cd4_t_cell"))

seur_obj <- AddMetaData(seur_obj, real_labels, col.name="real_labels")
pdf(file = "cd4_expr_real.pdf", width = 6.75, heigh = 4.2)
VlnPlot(seur_obj, features = c("CD4", "CD8A", "CD8B"), group.by = "real_labels", pt.size=0.0)
dev.off()

pdf(file = "cd4_expr_preds.pdf", width = 6.75, heigh = 4.2)
VlnPlot(seur_obj, features = c("CD4", "CD8A", "CD8B"), group.by = "sharp_preds", pt.size=0)
dev.off()


# "CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1"
pdf(file = "cd8_t_cell_expr.pdf", width = 8, heigh = 5.5)
VlnPlot(seur_obj, features = c("CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1"), group.by = "sharp_preds", pt.size=0)
dev.off()# "CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1"
