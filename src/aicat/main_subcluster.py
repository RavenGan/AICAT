import os
import json
from utils import calculate_diff_genes, sub_clustering, calculate_within_group_markers
import scanpy as sc
from AnnoCross import AnnCross_subcluster, refine_grouped_cluster_idx_dict


os.chdir("/Users/david/Desktop/CelltypeAgenticAI")

# set up the OpenAI API key CONFIDENTIAL
openai_api_key = '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'
species = "human"
tissue = 'Skin' #'Prostate' # 'Breast' #'Skin'
cell_type_col = 'Broad cell type'
res_save_path = "./res/2025_0507/skin"

# load the preprocessed data
data_path = f"./data/GTEx/processed/{tissue}_processed.h5ad"
adata = sc.read_h5ad(data_path)

### Perform differential expression analysis
# save_path_all_clusters = f"./data/GTEx/DE_tab/{tissue}_EachCluster_top_genes.json"
# cell_type_markers = calculate_diff_genes(adata, cell_type_col,
#                          save_path=save_path_all_clusters)

# ### In-depth annotation of single cell types===========================================
# anno_res = AnnoSingle(openai_api_key, cell_type_markers, species, tissue, res_save_path)
with open(f'{res_save_path}/AnnoSingle_skin_res_dict.json', 'r') as f:
    anno_res = json.load(f)

### Sub-clustering and annotation refinement===========================================
# Suppose we want to do sub-clustering and refine the annotation for a specific cluster
chosen_cluster_idx = "Fibroblast"
key_added = "subcluster"

# Perform sub-clustering for the chosen cluster
data_subcluster = sub_clustering(adata, 
                                 cell_type_col, 
                                 chosen_cluster_idx,
                                 key_added=key_added, # new cell type column name
                                 resolution=0.8)

# Save the subclustered data
save_path_subcluster = f"{res_save_path}/{tissue}_Subcluster_{chosen_cluster_idx}_top_genes.json"
subcluster_markers_dict = calculate_diff_genes(data_subcluster,
                                                key_added,
                                                save_path=save_path_subcluster)

subcluster_dict = AnnCross_subcluster(openai_api_key,
                                      anno_res,
                                      chosen_cluster_idx,
                                      subcluster_markers_dict,
                                      res_save_path,
                                      tissue,
                                      AnnoRound=0, 
                                      anno_level="celltype")

refine_group_dict = refine_grouped_cluster_idx_dict(subcluster_dict)

# exhaustive annotation refinement===========================================
AnnoRound = 1
while refine_group_dict != {}:
    refine_subcluster_markers_dict = calculate_within_group_markers(data_subcluster, 
                                                                      key_added,
                                                                      refine_group_dict,
                                                                      save_path=None)
    refine_subcluster_markers_dict = refine_subcluster_markers_dict[chosen_cluster_idx]
    subcluster_dict = AnnCross_subcluster(openai_api_key,
                                      anno_res,
                                      chosen_cluster_idx,
                                      refine_subcluster_markers_dict,
                                      res_save_path,
                                      tissue,
                                      AnnoRound=AnnoRound)
    AnnoRound += 1