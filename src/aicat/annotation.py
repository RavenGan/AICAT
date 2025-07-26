from .agent import CellTypeAgent, Annotate, group_clusters, cross_Annotate_single
from .flow import *
from .utils import chat_to_markdown

import os

def annotate_one_cluster(api_key, cluster_idx, marker_ls, species, 
                         tissue, res_save_path, initial_n_genes = 10, 
                         max_n_genes = 50, step_size = 10):
    """
    Annotate one cell cluster using the CellTypeAgent.

    Args:
        api_key (str): OpenAI API key.
        cluster_idx (str): Cluster index or name.
        marker_str (str): Comma-separated string of marker genes.
        species (str): Species name (e.g., "human").
        tissue (str): Tissue name (e.g., "prostate").
        res_save_path (str): Path to save the results.

    Returns:
        dict: Final cell type marker map.
    """
    # Set up the agent
    agent = CellTypeAgent(api_key=api_key, 
                          tools = [], #[wiki_tool],
                          ResponseFormat=Annotate,
                          verbose=False,
                          mode="single")

    # Check and create if the save path doesn't exist
    log_save_path = f"{res_save_path}/InDepthAnno/log_res"
    md_save_path = f"{res_save_path}/InDepthAnno/md_res"
    pdf_save_path = f"{res_save_path}/InDepthAnno/pdf_res"
    if not os.path.exists(log_save_path):
        os.makedirs(log_save_path)
    
    if not os.path.exists(md_save_path):
        os.makedirs(md_save_path)
    
    if not os.path.exists(pdf_save_path):
        os.makedirs(pdf_save_path)
    
    # Run the annotation flow
    annotation_track = annotation_flow(
        agent=agent,
        cluster_idx=cluster_idx,
        marker_ls=marker_ls,
        species=species,
        tissue=tissue,
        log_path=log_save_path,
        initial_n_genes=initial_n_genes,
        max_n_genes=max_n_genes,
        step_size=step_size

    )

    # Save the final cell type marker map
    print(f"The annotation tracking history is: {annotation_track}")

    # organize and save the chat history
    chat_hist = agent.chat_history
    chat_hist = chat_to_markdown(chat_hist)

    with open(f"{md_save_path}/{cluster_idx}_chat_history.md", "w") as f:
        f.write(chat_hist)
    # convert_md_to_pdf_with_pandoc(f"{res_save_path}/md_res/{cluster_idx}_chat_history.md", 
    #                             f"{res_save_path}/pdf_res/{cluster_idx}_chat_history.pdf")
    return annotation_track

def group_similar_clusters(api_key, anno_res, res_save_path,
                           anno_level = "celltype" # "lineage", "celltype", "fcn_anno"
                           ):

    """
    Group similar clusters based on their annotations.

    Args:
        api_key (str): OpenAI API key.
        anno_res (dict): Dictionary containing annotation results from previous results.

    Returns:
        dict: Dictionary with cluster indices grouped together.
    """
    # Set up the agent
    agent = CellTypeAgent(api_key=api_key, 
                          tools = [], #[wiki_tool],
                          ResponseFormat=group_clusters,
                          verbose=False,
                          mode="group_clusters")
    
    # Check and create if the save path doesn't exist
    log_save_path = f"{res_save_path}/GroupAnno/log_res"
    md_save_path = f"{res_save_path}/GroupAnno/md_res"
    pdf_save_path = f"{res_save_path}/GroupAnno/pdf_res"

    if not os.path.exists(log_save_path):
        os.makedirs(log_save_path)
        print(f"Created directory: {log_save_path}")
    
    if not os.path.exists(md_save_path):
        os.makedirs(md_save_path)
        print(f"Created directory: {md_save_path}")
    
    if not os.path.exists(pdf_save_path):
        os.makedirs(pdf_save_path)
        print(f"Created directory: {pdf_save_path}")
    
    # create cluster index-cluster id mapping
    # obtain list of dictionaries, keys are cluster ids, values are annotations.
    print("Organizing annotation results----")
    id_anno_dict_lst, cluster_id_dict = organize_anno_res(anno_res, 
                                                          anno_level=anno_level
                                                          )

    # obtained grouped cluster indices
    print("Grouping similar clusters----")
    grouped_cluster_idx_dict = group_anno(agent, id_anno_dict_lst, 
                                          cluster_id_dict, log_save_path)
    
    # organize and save the chat history
    chat_hist = agent.chat_history
    chat_hist = chat_to_markdown(chat_hist)

    print("Saving chat history----")
    with open(f"{md_save_path}/group_chat_history.md", "w") as f:
        f.write(chat_hist)

    return grouped_cluster_idx_dict


def annotation_cross(api_key, input_dict, 
                     res_save_path, group_id, AnnoRound):
    """
    Perform cross-annotation for a single cluster.

    Args:
        api_key (str): OpenAI API key.
        input_dict (dict): Dictionary containing input data for cross-annotation.
        res_save_path (str): Path to save the results.

    Returns:
        dict: Dictionary containing the cross-annotation results.
    """
    # Set up the agent
    agent = CellTypeAgent(api_key=api_key, 
                          tools = [], #[wiki_tool],
                          ResponseFormat=cross_Annotate_single,
                          verbose=False,
                          mode="cross")
    
    # Check and create if the save path doesn't exist
    log_save_path = f"{res_save_path}/CrossAnno/{group_id}_{AnnoRound}_log_res"
    md_save_path = f"{res_save_path}/CrossAnno/{group_id}_{AnnoRound}_md_res"
    pdf_save_path = f"{res_save_path}/CrossAnno/{group_id}_{AnnoRound}_pdf_res"
    
    if not os.path.exists(log_save_path):
        os.makedirs(log_save_path)
    
    if not os.path.exists(md_save_path):
        os.makedirs(md_save_path)
    
    if not os.path.exists(pdf_save_path):
        os.makedirs(pdf_save_path)

    # refine annotations for each cluster in a group
    print("----Refining annotations for each cluster in a group----")
    cross_single_dict = cross_anno_single(agent, input_dict, log_save_path)
    print("----Finished refining annotations for each cluster in a group----")

    # organize and save the chat history
    chat_hist = agent.chat_history
    chat_hist = chat_to_markdown(chat_hist)

    with open(f"{md_save_path}/{group_id}_cross_chat_history.md", "w") as f:
        f.write(chat_hist)
    
    return cross_single_dict


""" 
### unused function for cross-annotation and merging clusters

def annotation_cross_decide_merging(api_key, input_dict, res_save_path,
                                    decide_merging_anyway=False):

    agent = CellTypeAgent(api_key=api_key, 
                          tools = [], #[wiki_tool],
                          ResponseFormat=cross_Annotate_single,
                          verbose=False,
                          mode="cross")
    
    # Check and create if the save path doesn't exist
    log_save_path = f"{res_save_path}/CrossAnno/log_res"
    md_save_path = f"{res_save_path}/CrossAnno/md_res"
    pdf_save_path = f"{res_save_path}/CrossAnno/pdf_res"
    
    if not os.path.exists(log_save_path):
        os.makedirs(log_save_path)
    
    if not os.path.exists(md_save_path):
        os.makedirs(md_save_path)
    
    if not os.path.exists(pdf_save_path):
        os.makedirs(pdf_save_path)

    # refine annotations for each cluster in a group
    print("Refining annotations for each cluster in a group----")
    cross_single_dict = cross_anno_single(agent, input_dict)

    # separate clusters that need further check and those that are intact
    print("Separating clusters that need further check and those that are intact----")
    clusters_intact, clusters_further_check = organize_single_res(cross_single_dict,
                                                              input_dict)
    
    cluster_intact_output = {k: cross_single_dict[k] for k in clusters_intact if k in cross_single_dict}
    clusters_merged_output = {}

    if decide_merging_anyway:
        print("Decide to merge clusters anyway.-------")
        clusters_further_check = clusters_further_check + clusters_intact
        clusters_further_check_str = ", ".join(clusters_further_check)

        print(f"Check whether clusters {clusters_further_check_str} need to be merged----")

        # update the output format for cross annotation pairs
        agent.update_ResponseFormat(cross_Annotate_pair)

        # Decide whether to merge clusters or not
        cross_pair_dict = cross_anno_pair(agent, input_dict, clusters_further_check)

        pairs_checked = cross_pair_dict.keys()
        
        # keep all the decisions
        for pair in pairs_checked:
            clusters_merged_output[pair] = cross_pair_dict[pair]
    
    else:
        if len(clusters_further_check) >= 2:
            clusters_further_check_str = ", ".join(clusters_further_check)
            print(f"Check whether clusters {clusters_further_check_str} need to be merged----")
            
            # update the output format for cross annotation pairs
            agent.update_ResponseFormat(cross_Annotate_pair)

            # Decide whether to merge clusters or not
            cross_pair_dict = cross_anno_pair(agent, input_dict, clusters_further_check)

            pairs_checked = cross_pair_dict.keys()
            for pair in pairs_checked:
                # if the clusters are merged
                if cross_pair_dict[pair].merge_yes_or_no.lower() == "yes":
                    clusters_merged_output[pair] = cross_pair_dict[pair]
                
                # if the clusters are not merged
                elif cross_pair_dict[pair].merge_yes_or_no.lower() == "no":
                    target_idx, query_idx = pair
                    cluster_intact_output[target_idx] = cross_single_dict[target_idx]
                    cluster_intact_output[query_idx] = cross_single_dict[query_idx]
        else: 
            print("No clusters need to be merged.")
    
    # organize and save the chat history
    chat_hist = agent.chat_history
    chat_hist = chat_to_markdown(chat_hist)

    with open(f"{res_save_path}/md_res/cross_chat_history.md", "w") as f:
        f.write(chat_hist)

    return cluster_intact_output, clusters_merged_output
"""