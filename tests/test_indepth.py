from aicat.main_indepth import indepth_annotation
from dotenv import load_dotenv
import os

def test_annotate_celltypes():
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("Set OPENAI_API_KEY in your environment")

    adata_path = "tests/data/CRC_SMC05-T_processed.h5ad"
    species = "human"
    tissue = "primary colorectal cancer"
    cluster_col_name = 'Cell_type'
    data_name = "CRC_SMC05-T"
    save_path = f"tests/res/{data_name}" 

    # Call function
    indepth_annotation(api_key, 
                       adata_path,
                       species,
                       tissue,
                       cluster_col_name,
                       data_name=data_name,
                       save_path=save_path)

if __name__ == "__main__":
    test_annotate_celltypes()