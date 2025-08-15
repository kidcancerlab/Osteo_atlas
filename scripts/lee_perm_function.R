run_lee_perms <- function(sobj,
                          colname_1,
                          colname_2,
                          nsim = 10000,
                          nproc = 10,
                          k_nn = 5,
                          nb2listw_style = "C") {
    coords <- SeuratObject::GetTissueCoordinates(sobj)

    neighbors <-
        spdep::knearneigh(coords, k = k_nn, longlat = FALSE) |>
        spdep::knn2nb(row.names = rownames(coords)) |>
        spdep::nb2listw(style = nb2listw_style, zero.policy = FALSE)

    lee_test_res <- get_lee_p(
        nsim = nsim,
        cell_type_1 = sobj@meta.data[[colname_1]],
        cell_type_2 = sobj@meta.data[[colname_2]],
        neighbors = neighbors,
        nproc = nproc
    )

    return(lee_test_res)
}

get_lee_p <- function(nsim, cell_type_1, cell_type_2, neighbors, nproc) {
    real_lee <- spdep::lee(
        cell_type_1,
        cell_type_2,
        neighbors,
        n = length(cell_type_1)
    )$L

    sims <- parallel::mclapply(
        1:nsim,
        mc.cores = nproc,
        function(x) {
            # Permuting by shuffling the spatial location of the spots, but
            # keeping the data paired per
            # https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/how-bivariate-spatial-association-works.htm
            rand_order <- sample(seq_along(cell_type_1))
            spdep::lee(
                cell_type_1[rand_order],
                cell_type_2[rand_order],
                neighbors,
                n = length(cell_type_1)
            )$L
        }
    ) |>
        unlist()

    result_list <- get_perm_p(real_lee, sims)

    return(result_list)
}

get_perm_p <- function(real, sims) {
    nsim <- length(sims)

    # Negative lee values indicate a negative correlation
    # so if all the sims are greater than the real lee value, the p-value
    # should be less than the smallest possible p-value (1/nsim)
    if (real < 0) {
        if (all(sims > real)) {
            p_value <- paste0("p < ", 1 / nsim)
        } else {
            p_value <- sum(sims <= real) / nsim
        }
        # Positive lee values indicate a positive correlation
        # so if all the sims are less than the real lee value, the p-value
        # should be less than the smallest possible p-value (1/nsim)
    } else {
        if (all(sims < real)) {
            p_value <- paste0("p < ", 1 / nsim)
        } else {
            p_value <- sum(sims >= real) / nsim
        }
    }

    return(list(p_value = p_value, sims = sims, real_lee = real))
}