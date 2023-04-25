library(Seurat)

seur_obj <- readRDS("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/process_seur.rds")
seur_obj$sharp_preds[seur_obj$sharp_preds == 0] = "b_cells"
seur_obj$sharp_preds[seur_obj$sharp_preds == 1] = "cd14_monocytes"
seur_obj$sharp_preds[seur_obj$sharp_preds == 2] = "cd4_t_cell"
seur_obj$sharp_preds[seur_obj$sharp_preds == 3] = "cd56_nk"
seur_obj$sharp_preds[seur_obj$sharp_preds == 4] = "cd8_t_cell"
seur_obj$sharp_preds <- factor(seur_obj$sharp_preds, levels = c("b_cells", "cd14_monocytes", "cd56_nk", "cd8_t_cell", "cd4_t_cell"))

Idents(seur_obj) <- "sharp_preds"
print(FindMarkers(seur_obj, features = c("CD4"), ident.1 = "cd4_t_cell", min.pct = 0))

