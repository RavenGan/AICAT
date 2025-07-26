from .annotation import group_similar_clusters

import json

def AnnoGroup(openai_api_key, anno_res, tissue, res_save_path):
    """
    Group similar clusters based on annotation results.

    Args:
        openai_api_key (str): OpenAI API key.
        anno_res (dict): Annotation results from previous step.
        res_save_path (str): Path to save the results.

    Returns:
        dict: Dictionary of grouped clusters.
    """
    # Group similar clusters
    tissue = tissue.lower()
    
    res_dict = group_similar_clusters(openai_api_key, anno_res, res_save_path)

    with open(f"{res_save_path}/AnnoGroup_{tissue}_res_dict.json", "w") as f:
        json.dump(res_dict, f, indent=4)
    
    return res_dict
