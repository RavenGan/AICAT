from aicat.main_indepth import indepth_annotation
from dotenv import load_dotenv
import os

def test_annotate_celltypes():
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("Set OPENAI_API_KEY in your environment")

    adata_path = "tests/data/MCA_Spleen_processed.h5ad"
    species = "mouse"
    tissue = "spleen"
    cluster_col_name = 'Annotation'
    data_name = "MCA_Spleen"
    save_path = "tests/res" 

    # Call function
    indepth_annotation(api_key, 
                       adata_path,
                       species,
                       tissue,
                       cluster_col_name,
                       data_name=data_name,
                       save_path=save_path)
