from aicat.main_subcluster import subcluster_annotation
from dotenv import load_dotenv
import os

def test_subcluster_annotate():
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("Set OPENAI_API_KEY in your environment")

    adata_path = "tests/data/CRC_SMC05-T_processed.h5ad"
    tissue = "primary colorectal cancer"
    cluster_col_name = 'Cell_type'
    chosen_cluster = "B cells"
    AnnoSingle_res_path = "tests/res/CRC_SMC05-T/AnnoSingle_primary colorectal cancer_res_dict.json"
    save_path = f"tests/res/CRC_SMC05-T_subcluster" 

    # Call function
    subcluster_annotation(api_key, 
                          adata_path,
                          tissue,
                          cluster_col_name,
                          chosen_cluster,
                          AnnoSingle_res_path,
                          save_path=save_path,
                          key_added="subcluster",
                          resolution=0.8,
                          anno_level="celltype" # one of three levels of annotation
                          )

if __name__ == "__main__":
    test_subcluster_annotate()
