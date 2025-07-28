import argparse
import scanpy as sc
from .main_indepth import indepth_annotation

def main():
    parser = argparse.ArgumentParser(description="In-depth annotation of single cell types.")
    parser.add_argument("--openai_api_key", type=str, required=True, help="OpenAI API key for annotation.")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to the AnnData object.")
    parser.add_argument("--species", type=str, required=True, help="Species of the dataset (e.g., 'human', 'mouse').")
    parser.add_argument("--tissue", type=str, required=True, help="Tissue type of the dataset.")
    parser.add_argument("--cluster_col_name", type=str, required=True, help="Column name for cluster labels in AnnData.")
    parser.add_argument("--data_name", type=str, default=None, help="Optional name for the dataset.")
    parser.add_argument("--save_path", type=str, default=None, help="Directory to save results.")

    args = parser.parse_args()

    adata = sc.read(args.adata_path)

    indepth_annotation(args.openai_api_key,
                       adata,
                       args.species,
                       args.tissue,
                       args.cluster_col_name,
                       data_name=args.data_name,
                       save_path=args.save_path)