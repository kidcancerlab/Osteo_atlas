new_make_loom_files <- function(input_table,
         out_dir = "loom_output/",
         cluster_account,
         slurm_base = paste0(getwd(), "/slurmOut"),
         sbatch_base = "sbatch_") {
    
    rownames(input_table) <- input_table$sample_id

    #make out_dir end with a /
    out_dir <- ifelse(endsWith(out_dir, "/"),
                      out_dir,
                      paste0(out_dir, "/"))
    

    system(paste0("mkdir -p ", out_dir))

    #make temporary sbatch directory
    system(paste0("mkdir ", out_dir, "sbatch"))


    gtf_list <- list("mixed" = "10x-hg38-mm10",
                     "human" = "10x-hg38",
                     "mouse" = "10x-mm10")

    for (sid in rownames(input_table)) {
        #make sample output folders
        sid_out <- paste0(out_dir, sid)
        system(paste0("mkdir -p", sid_out))

        species <- input_table[sid, ]$species
        new_gene_path <- paste0(sid_out, "/genes.gtf.gz")
        system(paste0("cp /home/gdrobertslab/lab/GenRef/",
                      gtf_list[[species]],
                      "/genes/genes.gtf.gz ",
                      new_gene_path))

        #read in h5 object
        h5_object <- Read10X_h5(input_table[sid, ]$h5_path)
        #make sure only get gene expression data
        if (class(h5_object) == "list") {
            h5_object <- h5_object[["Gene Expression"]]
        }

        #make temporary directory with barcodes for current sample
        bcs <- colnames(h5_object)
        write.table(bcs,
                    paste0(sid_out, "/tmp_bcs.tsv"),
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE)

        #copy over bams
        old_bam_path <- input_table[sid, ]$bam_path
        tmp_bam_path <- paste0(sid_out, "/tmp_bam.bam")
        if (!file.exists(tmp_bam_path)) {
            system(paste("cp", old_bam_path, tmp_bam_path))
        }
    }
    #Make conda environment
    #going to make environment in location of R installation, it's likely
    #that people won't have any conda environments here
    #get location of rrrSingleCellUtils and append env name to it
    rrrscu <- find.package("rrrSingleCellUtils")
    env_path <- paste0(rrrscu, "/r_rna_velo")

    #check conda environment doesn't exist before creating
    conda_envs <- system("conda info --envs", intern = TRUE)

    if (sum(grepl(pattern = env_path, x = conda_envs)) == 0) {
        #only make conda environment if it doesn't already exist
        exists_conda <-
            system(paste0("conda env create -p ",
                          env_path,
                          " -f ",
                          paste0(rrrscu,
                                 "/make_environment.yml")))
    }

    #create bash array of sample_ids
    ids <- rownames(input_table)
    id_array <- paste(ids, collapse = " ")

    replace_tbl <-
        tibble::tribble(
            ~find,                      ~replace,
            "placeholder_account",      cluster_account,
            "placeholder_slurm_out",    paste0(out_dir, "/sbatch/out_"),
            "placeholder_slurm_error",  paste0(out_dir, "/sbatch/error_"),
            "placeholder_env_path",     env_path,
            "placeholder_max_array",    as.character(length(ids) - 1),
            "placeholder_id_array",     id_array,
            "placeholder_out_dir",      out_dir
        )

    use_sbatch_template(replace_tibble = replace_tbl,
                        template = paste0("../../rrrSingleCellUtils/inst/make_loom_files_2.sh"),
                        submit = TRUE,
                        file_dir = paste0(sbatch_dir, "/jobs"))
}

use_sbatch_template <- function(replace_tibble,
                                template,
                                file_dir = tempdir(),
                                temp_ext = ".sh",
                                temp_prefix = "sbatch_",
                                warning_label = "",
                                submit = TRUE) {
    sbatch_template <-
        readr::read_file(template)

    # Replace placeholders with real data
    for (i in seq_len(nrow(replace_tibble))) {
        sbatch_template <-
        stringr::str_replace_all(sbatch_template,
                                 pattern = replace_tibble$find[i],
                                 replacement = replace_tibble$replace[i])
    }

    temp_file <-
        tempfile(fileext = temp_ext,
                 tmpdir = file_dir,
                 pattern = temp_prefix)

    readr::write_file(sbatch_template, file = temp_file)

    if (submit == TRUE) {
        return_val <- system(paste("sbatch", temp_file))
    } else {
        return_val <- 0
    }

    if (return_val != 0) {
        stop(paste0(warning_label,
                    " sbatch submission failed. Error code ",
                    return_val))
    }
    return(0)
}
