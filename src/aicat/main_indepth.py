from .utils import calculate_diff_genes, calculate_within_group_markers
from .AnnoSingle import AnnoSingle
from .AnnoGroup import AnnoGroup
from .AnnoCross import AnnCross, refine_grouped_cluster_idx_dict

import os
import argparse
import scanpy as sc

def indepth_annotation(openai_api_key, 
                       adata_path,
                       species,
                       tissue,
                       cluster_col_name,
                       data_name=None,
                       save_path=None):
    """
    Main function to perform in-depth annotation of single cell types.
    This function orchestrates the entire annotation process including
    differential expression analysis, group clustering, and cross-group annotation.
    """
    try:
        if save_path is None:
            save_path = os.getcwd()
        elif os.path.exists(save_path):
            if not os.path.isdir(save_path):
                raise NotADirectoryError(f"{save_path} exists but is not a directory.")
        else:
            os.makedirs(save_path, exist_ok=True)
    except Exception as e:
        raise RuntimeError(f"Failed to set or create save path: {save_path}") from e
    
    # Load the AnnData object===========================================
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"AnnData file not found at {adata_path}")
    
    adata = sc.read_h5ad(adata_path)

    if cluster_col_name not in adata.obs.columns:
        raise KeyError(f"Column '{cluster_col_name}' not found in adata.obs. Available columns: {adata.obs.columns.tolist()}")

    # Perform differential expression analysis====================================
    if data_name is None:
        save_path_all_clusters = f"{save_path}/EachCluster_top_genes.json"
    else:
        save_path_all_clusters = f"{save_path}/{data_name}_EachCluster_top_genes.json"

    cell_type_markers = calculate_diff_genes(adata, 
                                             cluster_col_name,
                                             save_path=save_path_all_clusters)
    
    ### In-depth annotation of single cell types====================================
    anno_res = AnnoSingle(openai_api_key, 
                          cell_type_markers, 
                          species, 
                          tissue, 
                          save_path)
    
    ### Group similar clusters
    grouped_cluster_idx_dict = AnnoGroup(openai_api_key, 
                                         anno_res, 
                                         tissue, 
                                         save_path)
    
    # Perform differential expression analysis for cross-group======================
    if data_name is None:
        save_path_within_group = f"{save_path}/WithinGroup_top_genes.json"
    else:
        save_path_within_group = f"{save_path}/{data_name}_WithinGroup_top_genes.json"

    within_group_markers_dict = calculate_within_group_markers(adata, 
                                                               cluster_col_name,
                                                               grouped_cluster_idx_dict,
                                                               save_path=save_path_within_group)
    
    # Perform cross-group annotation===========================================
    cross_group_dict = AnnCross(openai_api_key, 
                                anno_res, 
                                within_group_markers_dict, 
                                grouped_cluster_idx_dict, 
                                save_path, 
                                tissue, 
                                AnnoRound=0)

    refine_group_dict = refine_grouped_cluster_idx_dict(cross_group_dict)

    # exhaustive annotation refinement===========================================
    AnnoRound = 1
    while refine_group_dict != {}:
        refine_within_group_markers_dict = calculate_within_group_markers(adata, 
                                                                          cluster_col_name,
                                                                          refine_group_dict,
                                                                          save_path=None)
        refine_cross_group_dict = AnnCross(openai_api_key, 
                                           anno_res, 
                                           refine_within_group_markers_dict, 
                                           refine_group_dict, 
                                           save_path, 
                                           tissue, 
                                           AnnoRound)
        refine_group_dict = refine_grouped_cluster_idx_dict(refine_cross_group_dict)
        AnnoRound += 1







# os.chdir("/Users/david/Desktop/CelltypeAgenticAI")
os.chdir("/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/CelltypeAgenticAI")

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--tissue", type=str, help="A chosen tissue from the data MCA")
args = parser.parse_args()

tissue = args.tissue

# set up the OpenAI API key CONFIDENTIAL
openai_api_key = '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'
species = "mouse"

cell_type_col = 'Annotation'
res_save_path = f"./res/2025_0708_MCA/{tissue}"

# Check if the directory exists, if not, create it
if not os.path.exists(res_save_path):
    os.makedirs(res_save_path)

# load the preprocessed data
data_path = f"./data/MCA_mouse/processed_SingleBatch/{tissue}_processed.h5ad"
adata = sc.read_h5ad(data_path)

