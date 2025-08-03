import argparse
from .main_scoring import scoring_annotation

def main():
    parser = argparse.ArgumentParser(description="Score annotations from AnnoSingle results. Require input organized from AnnoSingle_res_dict.json to CollectAnnoSingle_res_df.csv")
    parser.add_argument("--openai_api_key", type=str, required=True, help="OpenAI API key for scoring.")
    parser.add_argument("--tissue", type=str, required=True, help="Tissue type of the dataset.")
    parser.add_argument("--save_path_root", type=str, required=True, help="Root directory to save results.")
    parser.add_argument("--AnnoSingle_res_path", type=str, required=True, help="Path to the results of AnnoSingle.")

    args = parser.parse_args()

    scoring_annotation(args.openai_api_key,
                       args.tissue,
                       args.save_path_root,
                       args.AnnoSingle_res_path)
    
if __name__ == "__main__":
    main()