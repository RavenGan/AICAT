from .utils import calculate_diff_genes, calculate_within_group_markers
from .AnnoSingle import AnnoSingle
from .AnnoGroup import AnnoGroup
from .AnnoCross import AnnCross, refine_grouped_cluster_idx_dict

import os
import argparse
import scanpy as sc

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

# Perform differential expression analysis
save_path_all_clusters = f"{res_save_path}/{tissue}_EachCluster_top_genes.json"
cell_type_markers = calculate_diff_genes(adata, cell_type_col,
                         save_path=save_path_all_clusters)

### In-depth annotation of single cell types===========================================
anno_res = AnnoSingle(openai_api_key, cell_type_markers, species, 
                      tissue, res_save_path)
# with open(f'{res_save_path}/AnnoSingle_skin_res_dict.json', 'r') as f:
#     anno_res = json.load(f)

### Group similar clusters
grouped_cluster_idx_dict = AnnoGroup(openai_api_key, anno_res, tissue, 
                                     res_save_path)
# with open(f'{res_save_path}/AnnoGroup_skin_res_dict.json', 'r') as f:
#     grouped_cluster_idx_dict = json.load(f)

# Perform differential expression analysis for cross-group===========================
save_path_within_group = f"{res_save_path}/{tissue}_WithinGroup_top_genes.json"
within_group_markers_dict = calculate_within_group_markers(adata, 
                                                           cell_type_col,
                                                           grouped_cluster_idx_dict,
                                                           save_path=save_path_within_group)

# Perform cross-group annotation===========================================
cross_group_dict = AnnCross(openai_api_key, anno_res, 
                            within_group_markers_dict, 
                            grouped_cluster_idx_dict, 
                            res_save_path, tissue, AnnoRound=0)
# with open(f'{res_save_path}/skin/AnnoCross_skin_dict.json', 'r') as f:
#     cross_group_dict = json.load(f)

refine_group_dict = refine_grouped_cluster_idx_dict(cross_group_dict)

# exhaustive annotation refinement===========================================
AnnoRound = 1
while refine_group_dict != {}:
    refine_within_group_markers_dict = calculate_within_group_markers(adata, 
                                                            cell_type_col,
                                                            refine_group_dict,
                                                            save_path=None)
    refine_cross_group_dict = AnnCross(openai_api_key, anno_res, 
                                       refine_within_group_markers_dict, 
                                       refine_group_dict, 
                                       res_save_path, tissue, AnnoRound)
    refine_group_dict = refine_grouped_cluster_idx_dict(refine_cross_group_dict)
    AnnoRound += 1
