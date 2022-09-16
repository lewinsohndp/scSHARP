library(Seurat)

seur_obj <- readRDS("/home/groups/ConradLab/daniel/sharp_data/jung/process_seur.rds")
seur_obj$sharp_preds[seur_obj$sharp_preds == 0] = "spc"
seur_obj$sharp_preds[seur_obj$sharp_preds == 1] = "spd"
seur_obj$sharp_preds[seur_obj$sharp_preds == 2] = "spg"
seur_obj$sharp_preds <- factor(seur_obj$sharp_preds, levels = c("spg", "spc", "spd"))

# "mt-Rnr2", "Ncl", "Hsp90ab1", "Dazl", "Prrc2c"
pdf(file = "figures/spg_expr.pdf", width = 6.75, heigh = 4.2)
VlnPlot(seur_obj, features = c("mt.Rnr2", "Ncl","Hsp90ab1", "Dazl", "Prrc2c"), group.by = "sharp_preds", pt.size=0)
dev.off()

#"Ldhc", "Ubb", "Fabp9", "Pabpc1", "Meig1"
pdf(file = "figures/spc_expr.pdf", width = 6.75, heigh = 4.2)
VlnPlot(seur_obj, features = c("Ldhc", "Ubb", "Fabp9", "Pabpc1", "Meig1"), group.by = "sharp_preds", pt.size=0)
dev.off()

# "Tnp1", "Smcp", "Tsga8", "Dbil5", "Gm9999"
pdf(file = "figures/spd_expr.pdf", width = 6.75, heigh = 4.2)
VlnPlot(seur_obj, features = c("Tnp1", "Smcp", "Tsga8", "Dbil5", "Gm9999"), group.by = "sharp_preds", pt.size=0)
dev.off()

gene_set <- c("mt.Rnr2", "Ncl","Hsp90ab1", "Dazl", "Prrc2c", "Ldhc", "Ubb", "Fabp9", "Pabpc1", "Meig1", "Tnp1", "Smcp", "Tsga8", "Dbil5", "Gm9999")
expr = AverageExpression(seur_obj, features=gene_set, slot="scale.data", group.by="sharp_preds")
expr <- expr$RNA
expr <- expr[match(gene_set, row.names(expr)),]
write.csv(expr, "figures/testis_expr.csv")

seur_obj <- readRDS("/home/groups/ConradLab/daniel/sharp_data/pbmc_test/process_seur.rds")
seur_obj$sharp_preds[seur_obj$sharp_preds == 0] = "b_cells"
seur_obj$sharp_preds[seur_obj$sharp_preds == 1] = "cd14_monocytes"
seur_obj$sharp_preds[seur_obj$sharp_preds == 2] = "cd4_t_cell"
seur_obj$sharp_preds[seur_obj$sharp_preds == 3] = "cd56_nk"
seur_obj$sharp_preds[seur_obj$sharp_preds == 4] = "cd8_t_cell"
seur_obj$sharp_preds <- factor(seur_obj$sharp_preds, levels = c("b_cells", "cd14_monocytes", "cd56_nk", "cd8_t_cell", "cd4_t_cell"))

# "CD74", "HLA-DRA", "CD79A", "CD79B", "HLA-DRB1"
pdf(file = "figures/b_cell_expr.pdf", width = 8, height = 5.5)
VlnPlot(seur_obj, features = c("CD74", "HLA.DRA", "CD79A", "CD79B", "HLA.DRB1"), group.by = "sharp_preds", pt.size = 0)
dev.off()

# "FTH1", "TYROBP", "S100A9", "S100A8", "FTL"
pdf(file = "figures/cd14_monocyte_expr.pdf", width = 8, height = 5.5)
VlnPlot(seur_obj, features = c("FTH1", "TYROBP", "S100A9", "S100A8", "FTL"), group.by = "sharp_preds", pt.size = 0)
dev.off()

#"LTB", "FTH1", "CD3E", "MT-CO2","RPL34"
pdf(file = "figures/cd4_t_cell_expr.pdf", width = 8, height = 5.5)
VlnPlot(seur_obj, features = c("LTB", "FTH1", "CD3E", "MT.CO2","RPL34"), group.by = "sharp_preds", pt.size = 0)
dev.off()

# "GNLY", "NKG7", "TYROBP", "RPS21", "GZMB"
pdf(file = "figures/cd56_nk_expr.pdf", width = 8, height = 5.5)
VlnPlot(seur_obj, features = c("GNLY", "NKG7", "TYROBP", "RPS21", "GZMB"), group.by = "sharp_preds", pt.size = 0)
dev.off()

# "CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1"
pdf(file = "figures/cd8_t_cell_expr.pdf", width = 8, heigh = 5.5)
VlnPlot(seur_obj, features = c("CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1"), group.by = "sharp_preds", pt.size=0)
dev.off()

gene_set <- c("CD74", "HLA.DRA", "CD79A", "CD79B", "HLA.DRB1", "FTH1", "TYROBP", "S100A9", "S100A8", "FTL", "GNLY", "NKG7", "TYROBP", "RPS21", "GZMB", "CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1", "LTB", "FTH1", "CD3E", "MT.CO2","RPL34")
#gene_set <- c("CD74", "HLA.DRA", "CD79A", "CD79B", "HLA.DRB1", "FTH1", "TYROBP", "S100A9", "S100A8", "FTL", "LTB", "FTH1", "CD3E", "MT.CO2","RPL34", "GNLY", "NKG7", "TYROBP", "RPS21", "GZMB", "CD8B", "CD3E", "CTSW", "TPT1", "EEF1A1")
expr = AverageExpression(seur_obj, features=gene_set, slot="scale.data", group.by="sharp_preds")
expr <- expr$RNA
expr <- expr[match(gene_set, row.names(expr)),]
write.csv(expr, "figures/pbmc_expr.csv")
