library(rrrSingleCellUtils)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringr)
library(harmony)
library(cowplot)
library(ggsci)

set.seed(888)
setwd("~/../lab/Analysis/Ryan/RProjects/CellTypeAnnRefs/")

# Start the graphics device
png("Plots/HuOsteoPrimary%03d.png", width = 720, height = 720,
  type = "cairo")

# Start with GSE152048
# Make a list of sample names
s <- c("BC5", "BC6", "BC10", "BC11", "BC16",
       "BC17", "BC20", "BC21", "BC22")
qc <- c(18000, 25000, 25000, 30000, 70000,
        40000, 70000, 50000, 50000)
path <- c("Conventional", "Conventional", "Conventional",
          "Conventional", "Conventional", "Chondroblastic",
          "Chondroblastic", "Intraosseous", "Chondroblastic")
type <- c("Primary", "Primary", "Lung Met", "Primary", "Primary",
          "Lung Met", "Primary", "Primary", "Primary")

# Download, file, and extract the files from GEO
if(!dir.exists("PrimaryTumor/GSE152048")) {
  tar_dir <- "PrimaryTumor/GSE152048"
  dir.create(tar_dir)
  geo_pre <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/suppl/GSE152048_"
  for(i in seq_len(s)){
    gse_path <- str_c(geo_pre, s[i], ".matrix.tar.gz")
    tar_file <- str_c(tar_dir, "/", s[i], ".tar.gz ")
    download.file(gse_path, destfile = tar_file, method = "auto")
    untar(tar_file, exdir = tar_dir)
    file.remove(tar_file)
  }
}

# Create a vector that contains normalized Seurat objects for all GSE152048
#   samples
raw <- c()
p <- list()
for(i in 1:length(s)) {
  x <- tenx_load_qc(str_c("PrimaryTumor/GSE152048/",
                          s[i], "/"))
  x <- subset(x, subset = nCount_RNA < qc[i] & percent.mt <13)
  x$src <- s[i]
  x$type <- type[i]
  x$path <- path[i]
  x$gse <- "GSE152048"
  raw[[s[i]]] <- x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>%
    RunUMAP(dims = 1:20)
  p[[i]] <- (DimPlot(raw[[s[i]]], pt.size = 1, label = T, reduction = "umap") +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
}
rm (x)
plot_grid(plotlist = p, ncol = 3)

# Transition to download and process GSE162454
# Make a list of sample names
s <- c("OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6")
qc <- c(50000, 45000, 23000, 50000, 50000, 45000)
pre <- c("GSM4952363_", "GSM4952364_", "GSM4952365_", "GSM5155198_",
  "GSM5155199_", "GSM5155200_")

if(!dir.exists("PrimaryTumor/GSE162454")) {
  tar_dir <- "PrimaryTumor/GSE162454"
  dir.create(tar_dir)
  geo_pre <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162454/suppl/GSE162454_RAW.tar"
  tar_file <- str_c(tar_dir, "/", "GSE162454.tar.gz ")
  options(timeout = 300)
  download.file(geo_pre, destfile = tar_file, method = "auto")
  untar(tar_file, exdir = tar_dir)
  file.remove(tar_file)
  for(i in 1:length(s)) {
    bfile <- str_c(pre[i], s[i], "_barcodes.tsv.gz")
    ffile <- str_c(pre[i], s[i], "_features.tsv.gz")
    mfile <- str_c(pre[i], s[i], "_matrix.mtx.gz")
    samp_dir <- str_c(tar_dir, "/", s[i])
    dir.create(samp_dir)
    file.rename(str_c(tar_dir, "/", bfile),
      str_c(samp_dir, "/", "barcodes.tsv.gz"))
    file.rename(str_c(tar_dir, "/", ffile),
      str_c(samp_dir, "/", "features.tsv.gz"))
    file.rename(str_c(tar_dir, "/", mfile),
      str_c(samp_dir, "/", "matrix.mtx.gz"))
  }
}

# Create a vector that contains normalized Seurat objects for all GSE162454
#   samples
p <- list()
for(i in 1:length(s)) {
  x <- tenx_load_qc(str_c("PrimaryTumor/GSE162454/",
                          s[i], "/"))
  x <- subset(x, subset = nCount_RNA < qc[i] & percent.mt <18)
  x$src <- s[i]
  x$type <- "Primary"
  x$path <- "Conventional"
  x$gse <- "GSE162454"
  raw[[s[i]]] <- x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>%
    RunUMAP(dims = 1:20)
  p[[i]] <- (DimPlot(raw[[s[i]]], pt.size = 1, label = T, reduction = "umap") +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
}
rm(x)
plot_grid(plotlist = p, ncol = 3)

# Subsample and merge into one Seurat object and integrate with
# harmony (primary only)
prim <- c()
raw <- raw[!names(raw) %in% c("BC10", "BC17")]
for(i in 1:length(raw)) {
  prim[[i]] <- subset(raw[[i]], cells = sample(Cells(raw[[i]]), 1000))
}

comb <- merge(prim[[1]], y = prim[2:length(prim)],
                add.cell.ids = names(raw),
                project = "PrimaryReference")

comb <- comb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F)

comb <- RunHarmony(comb, group.by.vars = "src")
comb <- comb %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()

# Save a stopping point - merged and harmony aligned
qs::qsave(comb, "PrimaryTumor/comb.qs")
# Start from this stopping point
comb <- qs::qread("PrimaryTumor/comb.qs")

DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("Primary Tumors (integrated)") +
  theme(legend.position = "none")

DimPlot(comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("Primary Tumors (integrated)") 

DimPlot(comb, reduction = "umap", split.by = "src", ncol = 4) +
  coord_fixed() +
  ggtitle("Primary Tumors (integrated)")

## Optimizing UMAP parameters
# p <- list()
# for(i in 1:5) {
#   test <- RunUMAP(comb, reduction = "harmony", dims = 1:30,
#     set.op.mix.ratio = (i-1)*0.25)
#   test <- test %>%
#     FindNeighbors(reduction = "harmony", dims = 1:30) %>%
#     FindClusters()
#   p[[i]] <- DimPlot(test, reduction = "umap", label = T, repel = T) +
#       coord_fixed() +
#       ggtitle(str_c("set.op.mix.ratio = ", (i-1)*0.25)) +
#       theme(legend.position = "none")
# }
# plot_grid(plotlist = p, ncol = 3)

# p <- list()
# for(i in 1:5) {
#   test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
#     n.neighbors = (i * 10))
#   test <- test %>%
#     FindNeighbors(reduction = "harmony", dims = 1:30) %>%
#     FindClusters()
#   p[[i]] <- DimPlot(test, reduction = "umap", label = T, repel = T) +
#       coord_fixed() +
#       ggtitle(str_c("n.neighbors = ", (i * 10))) +
#       theme(legend.position = "none")
# }
# plot_grid(plotlist = p, ncol = 3)

# p <- list()
# for(i in 1:5) {
#   test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
#     n.epochs = (i * 100))
#   test <- test %>%
#     FindNeighbors(reduction = "harmony", dims = 1:30) %>%
#     FindClusters()
#   p[[i]] <- DimPlot(test, reduction = "umap", label = T, repel = T) +
#     coord_fixed() +
#     ggtitle(str_c("n.epochs = ", (i * 100))) +
#     theme(legend.position = "none")
# }
# plot_grid(plotlist = p, ncol = 3)

# p <- list()
# for(i in 1:5) {
#   test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
#     min.dist = (.00045 * 4 ^ i))
#   test <- test %>%
#     FindNeighbors(reduction = "harmony", dims = 1:30) %>%
#     FindClusters()
#   p[[i]] <- DimPlot(test, reduction = "umap", label = T, repel = T) +
#     coord_fixed() +
#     ggtitle(str_c("min.dist = ", (.00045 * 4 ^ i))) +
#     theme(legend.position = "none")
# }
# plot_grid(plotlist = p, ncol = 3)

# p <- list()
# for(i in 1:5) {
#   test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
#     local.connectivity = (2 ^ i))
#   test <- test %>%
#     FindNeighbors(reduction = "harmony", dims = 1:30) %>%
#     FindClusters()
#   p[[i]] <- DimPlot(test, reduction = "umap", label = T, repel = T) +
#     coord_fixed() +
#     ggtitle(str_c("local.connectivity = ", (2 ^ i))) +
#     theme(legend.position = "none")
# }
# plot_grid(plotlist = p, ncol = 3)

# Plot out cell cycle
comb <- kill_cc(comb)

# Begin process of identifying cell types and annotating the samples
# Calculate and show module scores from the paper
ms <- list()
ms$Osteoblastic <- c("RUNX2", "COL1A1", "CDH11", "IBSP")
ms$Chondroblastic <- c("SOX9", "ACAN", "PTH1R")
ms$Osteoclast <- c("ACP5", "CTSK", "MMP9")
ms$Myeloid <- c("CD74", "CD14", "FCGR3A")
ms$TCell <- c("CD3E", "IL7R", "CD8A", "CD4", "NKG7")
ms$NKCell <- c("NKG7", "GNLY")
ms$NKTCell <- c("NKG7", "GNLY", "CD3E")
ms$DCCell <- c("CD1C", "FCER1A", "CLEC9A", "CCR7", "CD14", "CD163")
ms$Fibroblast <- c("DCN", "COL1A1")
ms$Pericyte <- c("RGS5", "ACTA2")
ms$MSC <- c("MME", "THY1", "CXCL12", "SFRP2")
ms$Endothelial <- c("PECAM1", "VWF")
ms$Myoblast <- c("MYL1", "MYLPF")
ms$BCell <- c("MS4A1", "CD19", "JCHAIN")

mod_names <- names(ms)

g <- list()
for(i in 1:length(mod_names)) {
  comb <- AddModuleScore(comb, ms[i], name = names(ms[i]))
  g[[i]] <- FeaturePlot(comb, features = str_c(names(ms[i]), "1"), pt.size = 1,
    order = T, cols = c("lightgoldenrod", "darkred")) +
    coord_fixed()
}
plot_grid(plotlist = g, labels = LETTERS[1:length(g)], ncol = 4, nrow = 4)

# Perform automated cell type assignment
# Assign cell types (roughly, Human Primary Cell Atlas)
ref <- celldex::HumanPrimaryCellAtlasData()
comb_sce <- DietSeurat(comb) %>%
  as.SingleCellExperiment()

# Make cell type predictions using SingleR
comb_pred <- SingleR::SingleR(
  test = comb_sce,
  ref = ref,
  labels = ref$label.main)

# Transfer cell labels back to the Seurat object and plot
comb$pred1 <- comb_pred$labels
Idents(comb) <- comb$pred1
DimPlot(comb, pt.size = 1, label = T, repel = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Cell Type Predictions (Human Primary Cell Atlas)")
DimPlot(comb, pt.size = 1, label = T, repel = T, split.by = "src", ncol = 4) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Cell Type Predictions (Human Primary Cell Atlas) by sample")

# Assign cell types (roughly, Blueprint/ENCODE)
ref <- celldex::BlueprintEncodeData()

# Make cell type predictions using SingleR
comb_pred <- SingleR::SingleR(
  test = comb_sce,
  ref = ref,
  labels = ref$label.main)

# Transfer cell labels back to the Seurat object and plot
comb$pred2 <- comb_pred$labels
Idents(comb) <- comb$pred2
DimPlot(comb, pt.size = 1, label = T, repel = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Cell Type Predictions (Blueprint/ENCODE)")
DimPlot(comb, pt.size = 1, label = T, repel = T, split.by = "src", ncol = 4) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Cell Type Predictions (Blueprint/ENCODE) by sample")

# Save a stopping point - cell type predictions made
qs::qsave(comb, "PrimaryTumor/comb-pred.qs")
# Start from this stopping point
comb <- qs::qread("PrimaryTumor/comb-pred.qs")

# Generate gene marker lists for each cluster
Idents(comb) <- comb$seurat_clusters
p <- DimPlot(comb, pt.size = 1, label = T, repel = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Original Clusters")
print(p)

if (!file.exists("comb_marks.xlsx")){
  comb_marks <- FindAllMarkers(comb)
  xlsx::write.xlsx(comb_marks, "comb_marks.xlsx")
}

# Isolate the clusters that seem most consistent with tumor cells
if (!file.exists("selection.rds")) {
  selection <- CellSelector(p)
  saveRDS(selection, "selection.rds")
}
selection <- readRDS("selection.rds")
tumor <- subset(comb, cells=selection)

DimPlot(tumor, pt.size = 1, label = T, repel = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Isolated (likely) tumor cells")

# Regress out cell cycle for the tumor cells and replot
tumor_cc <- kill_cc(tumor, cc_regress = "Y")

DimPlot(tumor_cc, pt.size = 1, label = T, repel = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Tumor - CC regressed")

tumor_cc2 <- tumor_cc %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.6) %>%
    RunUMAP(dims = 1:20)

DimPlot(tumor_cc, pt.size = 1, label = T, repel = T) +
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("Tumor - CC regressed")

dev.off()
