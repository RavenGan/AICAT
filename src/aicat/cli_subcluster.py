import argparse
from .main_subcluster import subcluster_annotation

def main():
    parser = argparse.ArgumentParser(description="Sub-cluster annotation of single cell types.")
    parser.add_argument("--openai_api_key", type=str, required=True, help="OpenAI API key for annotation.")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to the AnnData object.")
    parser.add_argument("--tissue", type=str, required=True, help="Tissue type of the dataset.")
    parser.add_argument("--cluster_col_name", type=str, required=True, help="Column name for cluster labels in AnnData.")
    parser.add_argument("--chosen_cluster", type=str, required=True, help="Cluster to be sub-clustered.")
    parser.add_argument("--AnnoSingle_res_path", type=str, required=True, help="Path to the results of AnnoSingle.")
    parser.add_argument("--save_path", type=str, default=None, help="Directory to save results.")
    parser.add_argument("--key_added", type=str, default="subcluster", help="New cell type column name.")
    parser.add_argument("--resolution", type=float, default=0.8, help="Resolution for sub-clustering.")
    parser.add_argument("--anno_level", type=str, default="celltype", help="Annotation level.")

    args = parser.parse_args()

    subcluster_annotation(args.openai_api_key,
                          args.adata_path,
                          args.tissue,
                          args.cluster_col_name,
                          args.chosen_cluster,
                          args.AnnoSingle_res_path,
                          save_path=args.save_path,
                          key_added=args.key_added,
                          resolution=args.resolution,
                          anno_level=args.anno_level)
