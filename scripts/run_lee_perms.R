library(tidyverse)
library(Seurat)

source("scripts/lee_perm_function.R")

args <- commandArgs(trailingOnly = TRUE)
sample_num <- as.numeric(args[1])

print(sample_num)

cell_types <-
    qs::qread("output/spacexr/granular_references/cell_groups.qs") %>%
    unlist() %>%
    as.vector()

# We'll add the output to this table
cell_combos <-
    combn(cell_types, 2, simplify = TRUE) %>%
    t() %>%
    as.data.frame()

colnames(cell_combos) <- c("cell_type_1", "cell_type_2")

cell_combos$lee_stat <- NA
cell_combos$perm_p <- NA

spatial_sobj <-
    qs::qread("output/spatial_objects/spatial_list_level3_annotations.qs")[[sample_num]]

surrounding_cell_no <- 6
n_perms <- 100000

print(spatial_sobj)

for (i in seq_len(nrow(cell_combos))) {
    message(
        paste0(
            "Running Lee's Permutation Test for ",
            cell_combos$cell_type_1[i],
            " and ",
            cell_combos$cell_type_2[i],
            " in sample number ",
            sample_num
        )
    )

    set.seed(1337)

    lee_output <-
        run_lee_perms(
            spatial_sobj,
            colname_1 = cell_combos$cell_type_1[i],
            colname_2 = cell_combos$cell_type_2[i],
            nsim = n_perms,
            k_nn = surrounding_cell_no,
            nproc = 20
        )

    cell_combos$lee_stat[i] <- lee_output$real_lee
    cell_combos$perm_p[i] <- lee_output$p_value
}

write_tsv(
    cell_combos,
    paste0(
        "output/spacexr/granular_references/lee_perms/lee_perms_",
        spatial_sobj$sample_name[1],
        ".tsv"
    )
)

sessionInfo()
