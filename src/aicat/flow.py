from .steps import *
import logging

def setup_logger(loop_id, save_path):
    """
    Set up a logger for a loop with a unique ID.
    """
    # Create a new logger with a unique name
    logger = logging.getLogger(f"loop_logger_{loop_id}")
    logger.setLevel(logging.INFO)

    # Remove old handlers if they exist
    if logger.hasHandlers():
        logger.handlers.clear()

    # Console handler (for printing during runtime)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(console_handler)

    # File handler for this loop
    file_handler = logging.FileHandler(f"{save_path}/log_{loop_id}.log", mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    # Add handler to logger
    logger.addHandler(file_handler)

    return logger

"""
The prompt flow to annotate one cell cluster is as follows
|> step_lineage
    |> no: step_no_lineage ---> Abort/Ask
    |> yes: step_decide_one_lineage
        |> no: step_cannot_decide_one_lineage ---> Abort/Ask
        |> yes: step_celltypes
            |> no: step_no_celltypes ---> Abort/Ask
            |> yes: step_decide_one_celltype
                |> no: step_cannot_decide_one_celltype ---> Abort/Ask
                |> yes: step_fcn_anno
                    |> no: step_no_fcn_anno ---> Abort/Ask
                    |> yes: step_decide_one_fcn_anno
                        |> no: step_cannot_decide_one_fcn_anno
                        |> yes: ---> Final output
"""
def add_gene_until_success(agent, response_parsed, step_fn, marker_ls, n_gene_involved, 
                      step_size, max_n_genes, description, log, min_n_genes=10):
    """
    General loop for trying a step function until it returns 'yes' or the gene limit is reached.
    """
    n_gene_involved = max(n_gene_involved, min_n_genes)  # Ensure n_gene_involved is at least min_n_genes
    
    gene_ls = marker_ls[:n_gene_involved] # to avoid referencing issues
    gene_str = ", ".join(gene_ls)  # convert list to string
    
    while response_parsed.yes_or_no.lower() != "yes" and n_gene_involved < max_n_genes:
        log.info(f"-- Gene list cannot {description}. Add {step_size} more genes.")
        n_gene_involved += step_size
        if n_gene_involved > max_n_genes:
            n_gene_involved = max_n_genes
        gene_ls = marker_ls[:n_gene_involved]
        gene_str = ", ".join(gene_ls)
        log.info(f"Marker genes used: {gene_str}")
        response_parsed = step_fn(agent, response_parsed, gene_str)

    return response_parsed, gene_str, n_gene_involved


def annotation_flow(agent, cluster_idx, marker_ls, species, tissue, log_path,
                    initial_n_genes = 10, max_n_genes = 50, step_size = 10):
    
     # Start with the initial top genes. Will be updated in the flow
    n_gene_involved = initial_n_genes
    current_gene_ls = marker_ls[:n_gene_involved] 

    # set up the logger
    log = setup_logger(cluster_idx, log_path)
    log.info(f"Starting annotation flow for cluster {cluster_idx}.")
    log.info(f"Species: {species}, Tissue: {tissue}")
    log.info("-----------------------------------------")

    track_annotation_hierarchy = {}

    ### Step 1: Check if the lineage can be found-------
    # starting point of the flow
    log.info("Initial step: Checking lineage...")
    current_gene_str = ", ".join(current_gene_ls)  # convert list to string
    log.info(f"Marker genes used: {current_gene_str}")
    lineage_rep = step_lineage(agent, current_gene_str, species, tissue) 

    # iteratively check if the lineage is found and update the gene list if not
    lineage_rep, current_gene_str, n_gene_involved = add_gene_until_success(
        agent, lineage_rep, step_no_lineage, marker_ls, n_gene_involved, 
        step_size, max_n_genes, "identify cell lineage", log
    )
    
    # if the lineage still cannot be found, abort the flow
    if lineage_rep.yes_or_no.lower() != "yes":
        log.info("-- ❌ Cannot found any cell lineage given provided the gene list.")
        track_annotation_hierarchy["lineage"] = lineage_rep
        return track_annotation_hierarchy
    

    # if the lineages are found, proceed to decide one lineage
    log.info("-- ✅ Lineage found. Proceeding to decide one lineage...")
    log.info(f"Marker genes used: {current_gene_str}")
    decide_one_lineage = step_decide_one_lineage(agent, lineage_rep, current_gene_str)
    # iteratively check if one lineage can be found and update the gene list if not
    decide_one_lineage, current_gene_str, n_gene_involved = add_gene_until_success(
        agent, decide_one_lineage, step_decide_one_lineage, marker_ls, n_gene_involved, 
        step_size, max_n_genes, "decide one cell lineage", log
    )

    # if one single lineage still cannot be found, abort the flow
    if decide_one_lineage.yes_or_no.lower() != "yes":
        log.info("-- ❌ Cannot found one single cell lineage given provided the gene list.")
        track_annotation_hierarchy["lineage"] = decide_one_lineage
        return track_annotation_hierarchy
    
    # if one single lineage can be found, track and continue
    track_annotation_hierarchy["lineage"] = decide_one_lineage


    ### Step 2: Check if the cell type can be found-------
    log.info("-- ✅ One lineage determined. Proceeding to find cell types...")
    if n_gene_involved < 20:
        n_gene_involved = 20
    current_gene_ls = marker_ls[:n_gene_involved]
    current_gene_str = ", ".join(current_gene_ls)  # convert list to string
    log.info(f"Marker genes used: {current_gene_str}")
    celltype = step_celltypes(agent, decide_one_lineage, current_gene_str)

    # iteratively check if the cell types are found and update the gene list if not
    celltype, current_gene_str, n_gene_involved = add_gene_until_success(
        agent, celltype, step_celltypes, marker_ls, n_gene_involved, 
        step_size, max_n_genes, "identify cell types", log, min_n_genes=20 # ensure at least 20 genes are used for cell type annotation
    )
    # if the cell types still cannot be found, abort the flow
    if celltype.yes_or_no.lower() != "yes":
        log.info("-- ❌ Cannot found any cell types given provided the gene list and cell lineage.")
        track_annotation_hierarchy["celltype"] = celltype
        return track_annotation_hierarchy
    
    log.info("-- ✅ Cell types found. Proceeding to decide one single cell type...")
    log.info(f"Marker genes used: {current_gene_str}")
    decide_one_celltype = step_decide_one_celltype(agent, celltype, current_gene_str)
    # iteratively check if one cell type can be found and update the gene list if not
    decide_one_celltype, current_gene_str, n_gene_involved = add_gene_until_success(
        agent, decide_one_celltype, step_decide_one_celltype, marker_ls, n_gene_involved, 
        step_size, max_n_genes, "decide one cell type", log, min_n_genes=20
    )
    # if one single cell type still cannot be found, abort the flow
    if decide_one_celltype.yes_or_no.lower() != "yes":
        log.info("-- ❌ Cannot found one single cell type given provided the gene list and cell lineage.")
        track_annotation_hierarchy["celltype"] = decide_one_celltype
        return track_annotation_hierarchy
    
    # if one single cell type can be found, track and continue
    track_annotation_hierarchy["celltype"] = decide_one_celltype

    ### Step 3: Check if the functional annotations can be given-------
    log.info("-- ✅ One cell type determined. Proceeding to find functional annotations...")
    if n_gene_involved < 50:
        n_gene_involved = 50
    current_gene_ls = marker_ls[:n_gene_involved]
    current_gene_str = ", ".join(current_gene_ls)  # convert list to string
    log.info(f"Marker genes used: {current_gene_str}")
    fcn_annotations = step_fcn_anno(agent, decide_one_celltype, current_gene_str)

    # iteratively check if the functional annotations are found and update the gene list if not
    fcn_annotations, current_gene_str, n_gene_involved = add_gene_until_success(
        agent, fcn_annotations, step_fcn_anno, marker_ls, n_gene_involved, 
        step_size, max_n_genes, "find functional annotations", log, min_n_genes=50 # ensure at least 50 genes are used for cell type annotation
    )

    # if the functional annotations still cannot be found, abort the flow
    if fcn_annotations.yes_or_no.lower() != "yes":
        log.info("-- ❌ Cannot found any functional annotations given provided the gene list and cell lineage.")
        track_annotation_hierarchy["fcn_anno"] = fcn_annotations
        return track_annotation_hierarchy
    
    # if there are multiple functional annotations, decide one
    log.info("-- ✅ Functional annotations found. Proceeding to decide one annotation...")
    if len(fcn_annotations.cell_type) > 1:
        decide_one_fcn_anno = step_decide_one_fcn_anno(agent, fcn_annotations, current_gene_str)
        # iteratively check if one functional annotation can be found and update the gene list if not
        decide_one_fcn_anno, current_gene_str, n_gene_involved = add_gene_until_success(
            agent, decide_one_fcn_anno, step_decide_one_fcn_anno, marker_ls, n_gene_involved, 
            step_size, max_n_genes, "decide one functional annotation", log
        )
        # if one single functional annotation still cannot be found, abort the flow
        if decide_one_fcn_anno.yes_or_no.lower() != "yes":
            log.info("-- 🎉 There are multiple functional annotations given provided the gene list and cell lineage.")
            track_annotation_hierarchy["fcn_anno"] = decide_one_fcn_anno
            return track_annotation_hierarchy

        log.info("-- 🎉 One functional annotation found!")
        track_annotation_hierarchy["fcn_anno"] = decide_one_fcn_anno
        return track_annotation_hierarchy
    
    log.info("-- 🎉 One functional annotation found!")
    track_annotation_hierarchy["fcn_anno"] = fcn_annotations
    return track_annotation_hierarchy


"""
The prompt flow to group cell clusters with similar annotations is as follows
"""
# organize and extract the annotation results from previous flows
def organize_anno_res(anno_res, anno_level="celltype"):
    cluster_indices = anno_res.keys()
    cluster_id = 0

    # mapping a cluster index to a random cluster id
    # this is to avoid some cluster labels may indicate annotations.
    cluster_id_dict = {}
    
    # create the input for grouping the clusters
    id_anno_dict_lst = []
    
    for cluster_idx in cluster_indices:
        # get the annotation result for this cluster
        current_id = f"cluster_{cluster_id}"

        # in case the annotation result is empty or not a list
        try:
            anno_result = anno_res[cluster_idx][anno_level]["cell_type"][0]
        except (KeyError, IndexError, TypeError):
            anno_result = ""

        cluster_id_dict[cluster_idx] = current_id

        id_anno_dict_lst.append({
            "id": current_id,
            "cell_type": anno_result
        })

        cluster_id += 1
    
    return id_anno_dict_lst, cluster_id_dict

def group_anno(agent, id_anno_dict_lst, cluster_id_dict, log_path, retry_times=3):
    """
    Attempts to group annotated clusters using an LLM agent with retry logic.
    
    Parameters:
        agent: The LLM-based agent to use for annotation.
        id_anno_dict_lst: A list of annotation dictionaries.
        cluster_id_dict: A mapping from internal cluster IDs to annotation names.
        retry_times: Number of retry attempts in case of failure.
    
    Returns:
        grouped_cluster_idx_ls: A list of lists containing matched cluster indices per group.
    """
    grouped_cluster_idx_dict = {}
    response_parsed = None
    # set up the logger
    log_id = "group_anno"
    log = setup_logger(log_id, log_path)
    log.info(f"Group annotated clusters using an LLM agent with retry logic.")
    log.info("-----------------------------------------")

    for attempt in range(1, retry_times + 1):
        try:
            response_parsed = step_group_clusters(agent, id_anno_dict_lst)
            break  # Success, exit retry loop
        except Exception as e:
            log.info(f"[Retry {attempt}/{retry_times}] step_group_clusters failed: {e}")
            if attempt == retry_times:
                log.info("❌ All retry attempts failed. Returning empty result.")
                return grouped_cluster_idx_dict  # or raise an error if preferred

    # Parse and process the response
    try:
        dict_res = response_parsed.model_dump()
        results = dict_res.get("results", []) # Extract the results from the response

        group_num = 0
        for group in results:
            cluster_ids = group.get("group_ids", [])
            matched_cluster_idxs = [k for k, v in cluster_id_dict.items() if v in cluster_ids]
            grouped_cluster_idx_dict[f"group_{group_num}"] = matched_cluster_idxs
            group_num += 1
        log.info("Grouped clusters successfully.")

    except Exception as e:
        log.info(f"⚠️ Error while parsing or processing response: {e}")

    return grouped_cluster_idx_dict


"""
The prompt flow to perform cross cluster comparison is as follows
For each cell type in the pre-defined group:
    Run: step_cross_refine_single
    Output: Yes for refined cell type and its markers
            No for no refined cell type. Will further do cross comparison

For each pair of cell types with a response No:
    Run: step_cross_refine_pair
    Output: Yes for merged cell type and its markers
            No for no merged cell type. Keep the pair separate.

Make final decision on merged cell types and their markers.
"""
def create_subcluster_input_dict(anno_res, 
                                 chosen_cluster_idx,
                                 subcluster_markers_dict,
                                 anno_level="celltype"
                                 ):
    input_dict = {}
    sub_cluster_idx_ls = list(subcluster_markers_dict.keys())

    common_anno = anno_res.get(chosen_cluster_idx, {})
    common_anno_level_dict = common_anno.get(anno_level, {})
    cell_type_list = common_anno_level_dict.get("cell_type", [])
    anno_result = cell_type_list[0]

    common_marker_map = common_anno_level_dict.get("cell_type_marker_map", {})
    common_anno_markers = common_marker_map.get(anno_result, [])

    if not cell_type_list:
        print(f"⚠️ Warning: No annotation found for cluster {chosen_cluster_idx}. Skipping.")

    for sub_cluster_idx in sub_cluster_idx_ls:
        # Get the gene list for this cluster
        gene_list = subcluster_markers_dict[sub_cluster_idx]

        input_dict[sub_cluster_idx] = {
            "annotation": common_anno,
            "gene_list": ", ".join(gene_list),      # Convert list to string
            "marker_list": ", ".join(common_anno_markers)  # Convert list to string
        }
    
    return input_dict

def create_input_dict(anno_res, within_group_markers, 
                      grouped_cluster_idx_list,
                      anno_level="celltype"):
    """
    Create a dictionary for cross-annotation input.
    
    Parameters:
        anno_res: Dictionary containing annotation results. Generated from the annotation flow.
        within_group_markers: Dictionary containing gene markers for each cluster by within group comparisons
        grouped_cluster_idx_list: A list of grouped cluster indices.
        anno_level: The level of annotation to use (default is "celltype", other options are "lineage" and "fcn_anno").
    
    Returns:
        input_dict: A dictionary where keys are cluster indices and values are dictionaries with annotation details.
    """
    input_dict = {}

    for cluster_idx in grouped_cluster_idx_list:
        # Get the gene list for this cluster
        gene_list = within_group_markers.get(cluster_idx, [])

        # Get the annotation result for this cluster
        anno_cluster = anno_res.get(cluster_idx, {})
        anno_level_dict = anno_cluster.get(anno_level, {})
        cell_type_list = anno_level_dict.get("cell_type", [])

        if not cell_type_list:
            print(f"⚠️ Warning: No annotation found for cluster {cluster_idx}. Skipping.")
            continue

        anno_result = cell_type_list[0]

        # Get marker genes associated with the annotation
        marker_map = anno_level_dict.get("cell_type_marker_map", {})
        anno_markers = marker_map.get(anno_result, [])

        input_dict[cluster_idx] = {
            "annotation": anno_result,
            "gene_list": ", ".join(gene_list),      # Convert list to string
            "marker_list": ", ".join(anno_markers)  # Convert list to string
        }
    

    return input_dict


def cross_anno_single(agent, input_dict, log_path, retry_times=3):
    target_idx_ls = input_dict.keys()
    cross_single_dict = {}
    # set up the logger

    for target_idx in target_idx_ls:
        log = setup_logger(target_idx, log_path)
        log.info(f"--- Cross comparison: refinement for {target_idx} ---")
        log.info("-----------------------------------------")
        attempt = 0
        success = False

        while attempt < retry_times:
            try:
                response_parsed = step_cross_refine_single(agent, input_dict, target_idx)
                cross_single_dict[target_idx] = response_parsed
                log.info(f"✅ Attempt {attempt + 1} succeeded for {target_idx}.")
                success = True
                break  # exit retry loop on success
            except Exception as e:
                attempt += 1
                log.info(f"Attempt {attempt} failed for {target_idx}: {e}")

        if not success:
            log.info(f"❌ All {retry_times} attempts failed for {target_idx}. Skipping.")

    return cross_single_dict

def organize_single_res(cross_single_dict, input_dict):
    """
    Organize the results from cross_anno_single into a structured format.
    """
    target_idx_ls = input_dict.keys()
    # separate the clusters into two lists based on the response
    clusters_intact = []
    clusters_further_check = []

    for target_idx in target_idx_ls:
        cross_single_res = cross_single_dict[target_idx]
        if cross_single_res.yes_or_no.lower() == "yes":
            clusters_intact.append(target_idx)
        else:
            clusters_further_check.append(target_idx)

    return clusters_intact, clusters_further_check


def organize_group_res(cross_group_dict):
    """
    Organize group responses by extracting cluster IDs and their corresponding boolean decisions.

    Args:
        cross_group_dict (dict): Nested dictionary with group_id as keys and dictionaries of 
                                 cluster_id: {'yes_or_no': 'yes'/'no'} as values.

    Returns:
        dict: Dictionary mapping group_id to a tuple of (cluster_ids_list, yes_or_no_list),
              where yes_or_no_list contains booleans or None for unrecognized values.
    """
    organized_responses = {}
    for group_id, response_dict in cross_group_dict.items():
        cluster_ids = list(response_dict.keys())
        yes_or_no = [
            response.get('yes_or_no', '').strip().lower() == 'yes'
            if response.get('yes_or_no', '').strip().lower() in {'yes', 'no'} else None
            for response in response_dict.values()
        ]
        organized_responses[group_id] = (cluster_ids, yes_or_no)
    return organized_responses

def create_refine_dict(organized_responses):
    """
    Create a refined group dictionary by selecting cluster IDs 
    with a 'False' flag based on decision logic.

    Args:
        organized_responses (dict): Dictionary mapping group_id to a tuple 
                                    of (cluster_ids_list, yes_or_no_list).

    Returns:
        dict: Dictionary mapping group_id to a list of cluster_ids marked as False,
              but only if at least 2 are False and not all are False.
    """
    refine_group_dict = {}
    for group_id, (cluster_ids, decisions) in organized_responses.items():
        false_count = decisions.count(False)
        if 1 < false_count < len(decisions):
            sub_cluster_ids = [cid for cid, flag in zip(cluster_ids, decisions) if flag is False]
            refine_group_dict[group_id] = sub_cluster_ids

    return refine_group_dict
""" 
### unused function for cross-annotation and merging clusters

def cross_anno_pair(agent, input_dict, clusters_further_check, retry_times=3):
    
    # Perform pairwise cross comparison to decide merging clusters or not.
    # CHECK: len(clusters_further_check) >= 2

    pairs_to_check = list(itertools.combinations(clusters_further_check, 2))

    cross_pair_dict = {}

    for pair in pairs_to_check:
        target_idx, query_idx = pair
        print(f"--- Decide whether to merge: {target_idx} & {query_idx} ---")
        attempt = 0
        success = False

        while attempt < retry_times:
            try:
                response_parsed = step_cross_refine_pairwise(agent, input_dict, target_idx, query_idx)
                cross_pair_dict[(target_idx, query_idx)] = response_parsed
                print(f"✅ Attempt {attempt + 1} succeeded for {target_idx} & {query_idx}.")
                success = True
                break  # exit retry loop on success
            except Exception as e:
                attempt += 1
                print(f"Attempt {attempt} failed for {target_idx} & {query_idx}: {e}")

        if not success:
            print(f"❌ All {retry_times} attempts failed for {target_idx} & {query_idx}. Skipping.")
    return cross_pair_dict
""" 