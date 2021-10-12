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

# Define plot themes

gglayer_theme <- list(
  theme_minimal(),
  element_text(family = "DejaVu Sans"),
  plot.title = element_text(hjust = 0),
  coord_fixed()
)

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

# Create a vector that contains normalized Seurat objects for all GSE152048 samples
raw <- c()
for(i in 1:9) {
  x <- tenx_load_qc(str_c("PrimaryTumor/GSE152048/",
                          s[i], "/"))
  x <- subset(x, subset = nCount_RNA < qc[i] & percent.mt <13)
  x$src <- s[i]
  x$type <- type[i]
  x$path <- path[i]
  x$gse <- "GSE152048"
  x <- x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>%
    RunUMAP(dims = 1:20)
  print(DimPlot(x, pt.size = 1, label = T) +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
  raw <- c(raw, x)
}
rm(x)

# Transition to download and process GSE162454
# Make a list of sample names
s <- c("OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6")
qc <- c(50000, 45000, 23000, 50000, 50000, 45000)
pre <- c("GSM4952363_", "GSM4952364_", "GSM4952365_", "GSM5155198_", "GSM5155199_", "GSM5155200_")

if(!dir.exists("PrimaryTumor/GSE162454")) {
  tar_dir <- "PrimaryTumor/GSE162454"
  dir.create(tar_dir)
  geo_pre <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162454/suppl/GSE162454_RAW.tar"
  # gse_path <- str_c(geo_pre, s[i], ".matrix.tar.gz")
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
    file.rename(str_c(tar_dir, "/", bfile), str_c(samp_dir, "/", "barcodes.tsv.gz"))
    file.rename(str_c(tar_dir, "/", ffile), str_c(samp_dir, "/", "features.tsv.gz"))
    file.rename(str_c(tar_dir, "/", mfile), str_c(samp_dir, "/", "matrix.mtx.gz"))
  }
}

# Create a vector that contains normalized Seurat objects for all GSE162454 samples
for(i in 1:length(s)) {
  x <- tenx_load_qc(str_c("PrimaryTumor/GSE162454/",
                          s[i], "/"))
  x <- subset(x, subset = nCount_RNA < qc[i] & percent.mt <18)
  x$src <- s[i]
  x$type <- "Primary"
  x$path <- "Conventional"
  x$gse <- "GSE162454"
  x <- x %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>%
    RunUMAP(dims = 1:20)
  print(DimPlot(x, pt.size = 1, label = T) +
          coord_fixed() +
          theme(legend.position = "none") +
          ggtitle(str_c(s[i], " basic clustering")))
  raw <- c(raw, x)
}
rm(x)

# Subsample and merge into one Seurat object and integrate with harmony (primary only)
prim <- c()
for(i in 1:length(raw)) {
  x <- subset(raw[[i]], cells = sample(Cells(raw[[i]]), 1000))
  prim <- c(prim, x)
  rm(x)
}

comb <- merge(prim[[1]], y = prim[2:length(prim)],
                add.cell.ids = 1:length(prim),
                project = "PrimaryReference")

# Subset data to keep only primary tumors
comb <- subset(comb, type=="Primary")

comb <- comb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F)

comb <- RunHarmony(comb, group.by.vars = "src")
comb <- RunUMAP(comb, reduction = "harmony", dims = 1:30)
comb <- comb %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()

# Save a stopping point - merged and harmony aligned
qs::qsave(comb, "PrimaryTumor/comb.qs")
# Start from this stopping point
comb <- qs::qread("PrimaryTumor/comb.qs")

theme_set(theme_sf_light())

DimPlot(comb, reduction = "umap", label = T, repel = T) +
  coord_fixed() +
  ggtitle("Primary Tumors (integrated)") +
  theme(legend.position = "none") +
  theme_sf_light()

DimPlot(comb, reduction = "umap", group.by = "src") +
  coord_fixed() +
  ggtitle("Primary Tumors (integrated)") 

DimPlot(comb, reduction = "umap", split.by = "src", ncol = 4) +
  coord_fixed() +
  ggtitle("Primary Tumors (integrated)")

## Optimizing UMAP parameters
for(i in 1:5) {
  test <- RunUMAP(comb, reduction = "harmony", dims = 1:30,
    set.op.mix.ratio = (i-1)*0.25)
  test <- test %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters()
  p <- DimPlot(test, reduction = "umap", label = T, repel = T) +
    coord_fixed() +
    ggtitle(str_c("set.op.mix.ratio = ", (i-1)*0.25)) +
    theme(legend.position = "none")
  print(p)
}

for(i in 1:5) {
  test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
    n.neighbors = (i * 10))
  test <- test %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters()
  p <- DimPlot(test, reduction = "umap", label = T, repel = T) +
    coord_fixed() +
    ggtitle(str_c("n.neighbors = ", (i * 10))) +
    theme(legend.position = "none")
  print(p)
}

for(i in 1:5) {
  test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
    n.epochs = (i * 100))
  test <- test %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters()
  p <- DimPlot(test, reduction = "umap", label = T, repel = T) +
    coord_fixed() +
    ggtitle(str_c("n.epochs = ", (i * 100))) +
    theme(legend.position = "none")
  print(p)
}

for(i in 1:5) {
  test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
    min.dist = (.00045 * 4 ^ i))
  test <- test %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters()
  p <- DimPlot(test, reduction = "umap", label = T, repel = T) +
    coord_fixed() +
    ggtitle(str_c("min.dist = ", (.00045 * 4 ^ i))) +
    theme(legend.position = "none")
  print(p)
}

for(i in 1:5) {
  test <- RunUMAP(comb, reduction = "harmony", dims = 1:30, 
    local.connectivity = (2 ^ i))
  test <- test %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters()
  p <- DimPlot(test, reduction = "umap", label = T, repel = T) +
    coord_fixed() +
    ggtitle(str_c("local.connectivity = ", (2 ^ i))) +
    theme(legend.position = "none")
  print(p)
}

# # Ensure that groups are not just dominated by cell cycle effects
# comb <- kill_cc(comb)

# # Set a post-cc-regression stopping point
# qs::qsave(comb, "comb-ccReg.qs")
# # Start from this stopping point
# comb <- qs::qread("comb-ccReg.qs")

# DimPlot(primary, reduction = "umap", label = T, repel = T) +
#   coord_fixed() +
#   ggtitle("GSE152048 composite - primaries") +
#   theme(legend.position = "none")

# DimPlot(primary, reduction = "umap", group.by = "src") +
#   coord_fixed() +
#   ggtitle("GSE152048 composite - primaries") 

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
    coord_fixed() +
    theme(legend.position = "none")
}
plot_grid(plotlist = g, labels = LETTERS[1:length(g)], ncol = 4, nrow = 4)


# Perform automated cell type assignment
