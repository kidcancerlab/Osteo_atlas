library(tidyverse)
library(Seurat)

source("scripts/lee_perm_function.R")

args <- commandArgs(trailingOnly = TRUE)
combo_num <- as.numeric(args[1])

print(combo_num)

cell_types <-
    qs::qread("output/spacexr/granular_references/cell_groups.qs") %>%
    unlist() %>%
    as.vector()

cell_combos <-
    combn(cell_types, 2, simplify = TRUE) %>%
    t() %>%
    as.data.frame()

colnames(cell_combos) <- c("cell_type_1", "cell_type_2")

spatial_list <-
    qs::qread("output/spatial_objects/spatial_list_level3_annotations.qs")

surrounding_cell_no <- 6
n_perms <- 100000

results_table <-
    tibble(
        sample = names(spatial_list),
        cell_type_1 = cell_combos$cell_type_1[combo_num],
        cell_type_2 = cell_combos$cell_type_2[combo_num],
        lee_stat = NA,
        perm_p = NA
    )

for (i in seq_along(spatial_list)) {
    message(
        paste0(
            "Running Lee's Permutation Test for ",
            cell_combos$cell_type_1[combo_num],
            " and ",
            cell_combos$cell_type_2[combo_num],
            " in sample ",
            names(spatial_list)[i]
        )
    )

    set.seed(1337)

    lee_output <-
        run_lee_perms(
            spatial_list[[i]],
            colname_1 = cell_combos$cell_type_1[combo_num],
            colname_2 = cell_combos$cell_type_2[combo_num],
            nsim = n_perms,
            k_nn = surrounding_cell_no,
            nproc = 4
        )

    results_table$lee_stat[i] <- lee_output$real_lee
    results_table$perm_p[i] <- lee_output$p_value
}

write_tsv(
    results_table,
    paste0(
        "output/spacexr/granular_references/lee_perms/lee_perms_",
        combo_num,
        ".tsv"
    )
)
