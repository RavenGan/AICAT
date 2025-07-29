from .utils import calculate_diff_genes, sub_clustering, calculate_within_group_markers
from .AnnoCross import AnnCross_subcluster, refine_grouped_cluster_idx_dict
import os
import json
import scanpy as sc

def subcluster_annotation(openai_api_key, 
                          adata_path,
                          tissue,
                          cluster_col_name,
                          chosen_cluster,
                          AnnoSingle_res_path,
                          save_path=None,
                          key_added="subcluster",
                          resolution=0.8,
                          anno_level="celltype" # one of three levels of annotation
                          ):
    
    # Set up the save path if not provided or if it does not exist
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
    
    # Load the annotation results from AnnoSingle
    if not os.path.exists(AnnoSingle_res_path):
        raise FileNotFoundError(f"Annotation results file not found at {AnnoSingle_res_path}")

    with open(AnnoSingle_res_path, 'r') as f:
        anno_res = json.load(f)

    # Load the AnnData object===========================================
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"AnnData file not found at {adata_path}")
    
    adata = sc.read_h5ad(adata_path)

    if cluster_col_name not in adata.obs.columns:
        raise KeyError(f"Column '{cluster_col_name}' not found in adata.obs. Available columns: {adata.obs.columns.tolist()}")
    
    # Perform sub-clustering for the chosen cluster
    data_subcluster = sub_clustering(adata, 
                                     cluster_col_name, 
                                     chosen_cluster,
                                     key_added=key_added, # new cell type column name
                                     resolution=resolution)
    
    # Save the subclustered data
    save_path_subcluster = f"{save_path}/Subcluster_{chosen_cluster}_top_genes.json"
    subcluster_markers_dict = calculate_diff_genes(data_subcluster,
                                                   key_added,
                                                   save_path=save_path_subcluster)
    
    # Perform sub-cluster annotation
    subcluster_dict = AnnCross_subcluster(openai_api_key,
                                          anno_res,
                                          chosen_cluster,
                                          subcluster_markers_dict,
                                          save_path,
                                          tissue,
                                          AnnoRound=0, 
                                          anno_level=anno_level)
    
    refine_group_dict = refine_grouped_cluster_idx_dict(subcluster_dict)

    # exhaustive annotation refinement===========================================
    AnnoRound = 1
    while refine_group_dict != {}:
        refine_subcluster_markers_dict = calculate_within_group_markers(data_subcluster, 
                                                                        key_added,
                                                                        refine_group_dict,
                                                                        save_path=None)
        refine_subcluster_markers_dict = refine_subcluster_markers_dict[chosen_cluster]
        subcluster_dict = AnnCross_subcluster(openai_api_key,
                                              anno_res,
                                              chosen_cluster,
                                              refine_subcluster_markers_dict,
                                              save_path,
                                              tissue,
                                              AnnoRound=AnnoRound,
                                              anno_level=anno_level)
        AnnoRound += 1

