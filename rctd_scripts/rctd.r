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
        max_cores = 4,
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

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
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
