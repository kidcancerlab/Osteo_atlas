#!/usr/bin/env Rscript
library(rrrSingleCellUtils)
library(Seurat)
library(tidyverse)
library(spacexr)

# functions
run_rctd <- function(sp_ob, ref) {
    coords <- GetTissueCoordinates(sp_ob, image = "slice1") %>%
        dplyr::rename(x = imagerow, y = imagecol)
    #convert our object to an rctd object
    my_data <- spacexr::SpatialRNA(coords,
                                   GetAssayData(sp_ob, layer = "counts"))
    rctd_obj <- spacexr::create.RCTD(
        my_data,
        ref,
        max_cores = 1,
        UMI_min = 3,
        counts_MIN = 0,
        UMI_max = 900000000,
        CELL_MIN_INSTANCE = 0
    )
    rctd_out <- spacexr::run.RCTD(
        rctd_obj,
        doublet_mode = "doublet"
    )
    return(rctd_out)
}

size.factors <- list("OS1_Seurat" = 1500,
                     "OS2_Seurat" = 1400,
                     "OS3_Seurat" = 2700,
                     "OS4_Seurat" = 1700,
                     "OS5_Seurat" = 650,
                     "OS6_Seurat" = 3000,
                     "OS7_Seurat" = 2500,
                     "OS8_Seurat" = 1300)


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ob_name <- args[1]
ref_name <- args[2]

print(paste("Object is", ob_name))
print(paste("Reference is", ref_name))

# read in object list and select correct reference
spatial_list <- qs::qread("output/spatial_objects/unannotated_spatial_list.qs")
sp_ob <- spatial_list[[ob_name]]

# read in reference list and select correct reference
ref_list <- qs::qread("output/spacexr/granular_references/ref_list.qs")
ref <- ref_list[[ref_name]]

# Run rctd
rctd_res <- run_rctd(sp_ob = sp_ob, ref = ref)

# save results
qs::qsave(
    rctd_res,
    paste0(
        "output/spacexr/granular_references/",
        ref_name,
        "/",
        ob_name,
        ".qs"
    )
)

# add results to object for plotting
norm_weights <- spacexr::normalize_weights(rctd_res@results$weight)
sp_ob <- AddMetaData(sp_ob, norm_weights)

# get cell types from reference object
cell_types <- unique(ref@cell_types)

# get matrix of deconvolution scores
ann_mat <- sp_ob@meta.data[ , cell_types]

# plot heatmap of deconvolution score correlation between all cell types
pdf(
    paste0(
        "output/figures/spatial/spacexr/granular_",
        ref_name,
        "/",
        ob_name,
        "_heatmap.pdf"
    )
    width = 900,
    height = 900
)
pheatmap::pheatmap(
    cor(ann_mat, method = "spearman"),
    scale = "none",
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    angle_col = 45,
    fontsize = 16)
dev.off()

# make pdf of all rctd scores w h&e and cumulative tumor score at the start
sp_ob$Tumor_Cumulative <-
    sp_ob$Basal_Progenitor +
    sp_ob$Fibrogenic +
    sp_ob$Interactive +
    sp_ob$MP_Progenitor +
    sp_ob$Proliferative +
    sp_ob$Stressed
pdf(
    paste0("output/figures/spatial/spacexr/granular_",
    ref_name,
    "/",
    ob_name,
    "_scores.pdf")
    height = 4,
    width = 4
)
# H&E first
print(
    SpatialDimPlot(
        sp_ob,
        pt.size.factor = 0) + NoLegend()
)

# Now Tumor Cumulative
print(
    SpatialFeaturePlot(
        sp_ob,
        features = "Tumor_Cumulative",
        pt.size.factor = size.factors[[ob_name]],
        image.alpha = 0
    )
)

# now all cell type scores
for (ct in cell_types) {
    print(
        SpatialFeaturePlot(
            sp_ob,
            features = ct,
            pt.size.factor = size.factors[[ob_name]],
            image.alpha = 0
        )
    )
}
dev.off()