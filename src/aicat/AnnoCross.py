from .flow import create_subcluster_input_dict, create_input_dict, organize_group_res, create_refine_dict
from .annotation import annotation_cross

import json

def AnnCross_subcluster(api_key,
                        anno_res,
                        chosen_cluster_idx,
                        subcluster_markers_dict,
                        res_save_path,
                        tissue,
                        AnnoRound=0, 
                        anno_level="celltype"):
    input_dict = create_subcluster_input_dict(anno_res, 
                                              chosen_cluster_idx,
                                              subcluster_markers_dict,
                                              anno_level=anno_level)
    
    # Perform sub-cluster comparison and decide whether to merge sub-clusters
    subcluster_dict = annotation_cross(api_key, input_dict, res_save_path,
                                       group_id=chosen_cluster_idx,
                                       AnnoRound=AnnoRound)
    
    # convert Annotate object to dict for serialization
    subcluster_dict = {k: v.model_dump() for k, v in subcluster_dict.items()}

    # make the format the same as AnnCross
    subcluster_dict = {chosen_cluster_idx: subcluster_dict}

    with open(f"{res_save_path}/AnnCrossSubCluster_{tissue}_{chosen_cluster_idx}_round{AnnoRound}_dict.json", "w") as f:
        json.dump(subcluster_dict, f, indent=4)

    return subcluster_dict
    
    

def AnnCross_OneGroup(api_key, 
                      anno_res,
                      group_id,
                      within_group_markers, 
                      grouped_cluster_idx_list,
                      res_save_path,
                      AnnoRound,
                      anno_level="celltype"):
    
    #  prepare the input dictionary for cross-cluster comparison    
    input_dict = create_input_dict(anno_res, 
                                   within_group_markers, 
                                   grouped_cluster_idx_list,
                                   anno_level=anno_level)
    
    # Perform cross-cluster comparison and decide whether to merge clusters
    cross_single_dict = annotation_cross(api_key, input_dict, 
                                         res_save_path, group_id, AnnoRound)
    
    # convert Annotate object to dict for serialization
    cross_single_dict = {k: v.model_dump() for k, v in cross_single_dict.items()}

    return cross_single_dict

def AnnCross(api_key, anno_res, within_group_markers_dict, grouped_cluster_idx_dict, 
             res_save_path, tissue, AnnoRound=0, anno_level="celltype"):
    """
    Perform cross-cluster annotation and save the results.

    Args:
        api_key (str): OpenAI API key.
        anno_res (dict): Annotation results from previous step.
        within_group_markers_dict (dict): Dictionary of within-group markers.
        grouped_cluster_idx_dict (dict): Dictionary of grouped cluster indices.
        res_save_path (str): Path to save the results.
        tissue (str): Tissue name.

    Returns:
        dict: Dictionary of cross-cluster annotation results.
    """
    tissue = tissue.lower()
    cross_group_dict = {}

    # Obtain all the group IDs
    group_ids = list(within_group_markers_dict.keys())

    for group_id in group_ids:
        print(f"Processing group: {group_id}")
        within_group_markers = within_group_markers_dict[group_id]
        grouped_cluster_idx_list = grouped_cluster_idx_dict[group_id]

        # Perform cross-cluster comparison and decide whether to merge clusters
        cross_single_dict = AnnCross_OneGroup(api_key, 
                                              anno_res,
                                              group_id,
                                              within_group_markers, 
                                              grouped_cluster_idx_list,
                                              res_save_path,
                                              AnnoRound,
                                              anno_level=anno_level)
        
        cross_group_dict[group_id] = cross_single_dict
    
    with open(f"{res_save_path}/AnnoCross_{tissue}_round{AnnoRound}_dict.json", "w") as f:
        json.dump(cross_group_dict, f, indent=4)
    
    return cross_group_dict

def refine_grouped_cluster_idx_dict(cross_group_dict):
    organize_responses = organize_group_res(cross_group_dict)
    refine_group_dict = create_refine_dict(organize_responses)
    return refine_group_dict


