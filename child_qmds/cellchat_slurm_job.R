## Celchat with  shared genes
library(monocle3)
library(rrrSingleCellUtils)
library(rrrSnvs)
library(Seurat)
library(ggrepel)
library(tidyverse)
library(harmony)
library(cowplot)
library(clustree)
library(data.table)
library(hdf5r)
library(Rmagic)
library(knitr)
library(scATOMIC)
library(copykat)
library(SCENIC)
library(stashPlot)
library(ggalluvial)
library(RColorBrewer)
library(gridExtra)
library(CellChat)
library(patchwork)
library(circlize)
library(NMF)

# Pull command-line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_number <- args[1]
group <- args[2]

# function to run cellchat
run_cellchat <- function(sobject,
                        group.by="celltype",
                        species_db="human") {
    data.input <- GetAssayData(sobject, assay="RNA", slot="data")  # log-normalized
    meta <- data.frame(labels = sobject[[group.by]][,1], row.names=colnames(sobject))
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    cellchat@DB <- if (species_db=="human") CellChatDB.human else CellChatDB.mouse
    future::plan("multisession", workers = 10) # do parallel
    cellchat <- 
        subsetData(cellchat) %>%
        identifyOverExpressedGenes() %>%
        identifyOverExpressedInteractions() %>%
        computeCommunProb() %>%
        filterCommunication(min.cells = 10) %>%
        computeCommunProbPathway() %>%
        aggregateNet() %>%
        netAnalysis_computeCentrality()
        # computeNetSimilarity(type = "functional") %>%
        # netEmbedding(type = "functional") %>%
        # netClustering(type = "functional", do.parallel = FALSE) %>%
        # computeNetSimilarity(type = "structural") %>%
        # netEmbedding(type = "structural") %>%
        # netClustering(type = "structural", do.parallel = FALSE)
  return(cellchat)
}

# load the shared genes
shared_genes <-
    read_tsv("input/shared_genes_across_groups_for_cellchat.txt",
            col_names = "gene")  %>%
            as.data.frame()

object <-
    qs::qread(str_c("output/seurat_objects/final_tumor_vs_stroma/",
                    group,
                    ".qs"))
object <-
    subset(object, features = shared_genes$gene) %>%
    process_seurat()

samples <- unique(object$sample_name)
sample <-  samples[as.numeric(sample_number)]

warning(str_c("Processing sample ", sample, " of group ", group,
            " which is sample number ", sample_number, " of ", length(samples), " samples in this group"))

if (!file.exists(str_c("output/cellchat_objects_sharedgenes/", group, "/", sample, "_AnnL2_whole", ".qs"))) {
    sub_object <- subset(object, sample_name == sample)
    # Subset to max 500 cells per celltype (Ann_Level2)
    celltype_counts <- table(sub_object$Ann_Level3)
    sub_object <- SetIdent(sub_object, value = "Ann_Level3")
    cells_to_keep <- unlist(lapply(names(celltype_counts), function(ct) {
        ct_cells <- WhichCells(sub_object, idents = ct)
        if (length(ct_cells) > 100) {
            sample(ct_cells, 100)
        } else {
            ct_cells
        }
    }))
    sub_object <- 
        subset(object, cells = cells_to_keep)
    
    cellchat <-
        run_cellchat(sobject = sub_object,
                    group.by="Ann_Level2",
                    species_db=sub_object$organism[1])
    if (!dir.exists(str_c("output/cellchat_objects_sharedgenes", "/", group))) {
        dir.create(str_c("output/cellchat_objects_sharedgenes", "/", group),
                    recursive = TRUE)
    }
    qs::qsave(cellchat, str_c("output/cellchat_objects_sharedgenes/", group, "/", sample, "_AnnL2_whole", ".qs"))
} else {
    print(str_c("output/cellchat_objects_sharedgenes/", group, "/", sample, "_AnnL2_whole", ".qs already exists"))
}
