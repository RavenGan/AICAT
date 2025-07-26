from .annotation import annotate_one_cluster

import json

def AnnoSingle(api_key, cell_type_markers, species, tissue, res_save_path):
    """
    Annotate single cell types using the provided markers.

    Args:
        api_key (str): OpenAI API key.
        cell_type_markers (dict): Dictionary of cell type markers.
        species (str): Species name.
        tissue (str): Tissue name.
        res_save_path (str): Path to save the results.

    Returns:
        None
    """
    species = species.lower()
    tissue = tissue.lower()
    res_dict = {}

    for cluster_idx, marker_ls in cell_type_markers.items():
        cluster_idx = cluster_idx.replace("/", "_")  # replace "/" with "_" to avoid file path issues
        print("========================================")
        print(f"Cluster: {cluster_idx}")

        for attempt in range(3):  # try three times if there's an error
            try:
                res = annotate_one_cluster(api_key=api_key,
                                           cluster_idx=cluster_idx,
                                           marker_ls=marker_ls,
                                           species=species,
                                           tissue=tissue,
                                           res_save_path=res_save_path,
                                           initial_n_genes=10,
                                           max_n_genes=60,
                                           step_size=10)

                # convert Annotate object to dict for serialization
                res_serialized = {k: v.model_dump() for k, v in res.items()}
                res_dict[cluster_idx] = res_serialized
                print(f"Cluster {cluster_idx} annotated successfully.")
                break  # success, exit retry loop

            except Exception as e:
                print(f"Attempt {attempt + 1} failed for cluster {cluster_idx}: {e}")
                if attempt == 1:
                    print(f"Error annotating cluster {cluster_idx} after retry.")

    with open(f"{res_save_path}/AnnoSingle_{tissue}_res_dict.json", "w") as f:
        json.dump(res_dict, f, indent=4)
    
    return res_dict
