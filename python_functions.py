import os
import re
import anndata
import pandas as pd
import scvelo as scv
import numpy as np
import matplotlib as plt
plt.use("Agg")
from matplotlib import pyplot
import scanpy as py


### Python Functions

# `loom_to_an` will take all loom files associated with a seurat object, load them in one-by-one, add in appropriate metadata, and then merge them into a single anndata object
# `calc_velo` will calculate the rna velocity information from our loom files
# `write_obs` will take the output of loom_to_an, pass it to calc_velo, and then write it off
# `save_plots` will make stream plots for fdl and umap and group cells by seurat cluster and annotation, and write off each plot into output/figures/rna_velocity/


def loom_to_an(obj_name, loom_dir, metadata_dir):
    #Get sample id's from metadata files
    md_files=os.listdir(metadata_dir + "/" + obj_name)
    #initialize list to hold anndata objects
    ad_obs = []
    for one_md in md_files:
        #get sample id from file name
        sample_id = re.sub("_metadata.csv", "", one_md)
        #add object to list
        tmp = anndata.read_loom(loom_dir + "/" + sample_id + "/" + sample_id + ".loom")
        #read in metadata file
        md_path = metadata_dir + "/" + obj_name + "/" + one_md
        metadata = pd.read_csv(md_path, index_col = "bc")
        #Subset anndata to only barcodes in metadata folder
        tmp = tmp[metadata.index]
        tmp.obs = metadata.filter(regex = '^(?!umap_|PC_|fdl_|harmony_).*$').reindex(tmp.obs.index)
        #Add pca embeddings
        tmp.obsm["X_pca"] = metadata.filter(regex = "PC_").reindex(tmp.obs.index).to_numpy()
        #Add umap embeddings
        tmp.obsm["X_umap"] = metadata.filter(regex = "umap_").reindex(tmp.obs.index).to_numpy()
        #add fdl embeddings
        tmp.obsm["X_fdl"] = metadata.filter(regex = "fdl_").reindex(tmp.obs.index).to_numpy()
        #add harmony embeddings
        tmp.obsm["X_harmony"] = metadata.filter(regex = "harmony_").reindex(tmp.obs.index).to_numpy ()
        #make var names unique (this step is necessary to merge anndata list)
        tmp.var_names_make_unique()
        #add object to list
        ad_obs.append(tmp)
    #merge anndata list
    merged_anndata = anndata.concat(ad_obs)
    #return merged anndata
    return merged_anndata


def calc_velo(ad_ob, mode = "dynamical"):
    scv.pp.filter_and_normalize(ad_ob, n_top_genes = 2000)
    sc.pp.neighbors(ad_ob, n_pcs = 30, n_neighbors = 30)
    scv.pp.moments(ad_ob)
    if mode == "dynamical":
        scv.tl.recover_dynamics(ad_ob, n_jobs = 30)
    scv.tl.velocity(ad_ob, mode = mode)
    scv.tl.velocity_graph(ad_ob, backend = "threading", n_jobs = 30)

def write_obs(ob, method, mode = "dynamical", out_prefix = "loom_output/split_ad/"):
    # make out_prefix end with a "/"
    if not out_prefix.endswith('/'):
        out_prefix = out_prefix + "/"
    print("starting on " + ob)
    merged_ad=loom_to_an(obj_name = ob,
                     loom_dir = "loom_output/samples",
                     metadata_dir = "loom_output/split_metadata/" + method)
    print("made " + ob + " , now calculating velocity")
    calc_velo(merged_ad, mode = mode)
    if not os.path.isdir(out_prefix + method):
        os.makedirs(out_prefix + method)
    out_path = out_prefix + method + "/" + ob + ".ad"
    merged_ad.write(filename = out_path)
    print("calculated velocity and saved off " + ob + " for " + method)


def save_plots(ob, method, file_prefix = "output/figures/rna_velocity/", ad_prefix = "loom_output/split_ad/sc/"):
    merged_ad = anndata.read(ad_prefix + ob + ".ad")
    file_prefix = file_prefix + method + "/"
    # make file prefix end with "/"
    if not file_prefix.endswith('/'):
        out_prefix = out_prefix + "/"
    # check file prefix directories exist, and if not create them
    if not os.path.isdir(file_prefix):
        os.makedirs(file_prefix)
    #remove "Tumor_" from ann_level 3
    merged_ad.obs["Ann_Level3"] = merged_ad.obs["Ann_Level3"].str.replace("Tumor_", "")
    #Now start making plots
    scv.pl.velocity_embedding_stream(merged_ad,
                                     basis = "umap",
                                     color = "seurat_clusters",
                                     show = False,
                                     palette = hex_codes,
                                     sort_order = False,
                                     figsize = (7, 7),
                                     legend_loc = "none")
    pyplot.savefig(file_prefix + ob + "_umap_clusters.svg")
    #seurat_clusters in fdl
    scv.pl.velocity_embedding_stream(merged_ad,
                                     basis = "fdl",
                                     color = "seurat_clusters",
                                     palette = hex_codes,
                                     show = False,
                                     sort_order = False,
                                     figsize = (7, 7),
                                     legend_loc = "none")
    pyplot.savefig(file_prefix + ob + "_fdl_clusters.svg")
    #Now plot cell types
    scv.pl.velocity_embedding_stream(merged_ad,
                                     basis = "umap",
                                     color = "Ann_Level3",
                                     palette = color_mapping,
                                     show = False,
                                     sort_order = False,
                                     figsize = (7, 7),
                                     legend_loc = "none")
    pyplot.savefig(file_prefix + ob + "_umap_annotations.svg")
    #cell types in fdl
    scv.pl.velocity_embedding_stream(merged_ad,
                                     basis = "fdl",
                                     color = "Ann_Level3",
                                     palette = color_mapping,
                                     show = False,
                                     sort_order = False,
                                     figsize = (7, 7),
                                     legend_loc = "none")
    pyplot.savefig(file_prefix + ob + "_fdl_annotations.svg")
