### ===========================================================================
### steps used for single cluster annotations
def step_lineage(agent, marker_str, species, tissue):
    """
    Step: Initial query to the agent to identify broad cell lineages based on marker genes.
    Args:
        agent: The AI agent for annotation.
        marker_str (str): Comma-separated string of marker genes.
        species (str): Species name (e.g., "human").
        tissue (str): Tissue name (e.g., "prostate"). 
    """
    query = (
        "This is a list of top differential genes for a cell cluster "
        f"from {species} {tissue} tissues. The genes are: {marker_str}. "
        "Does the list contain lineage markers? Focus on decisive marker genes for cell lineage identification. "
        "If yes, decide the major cell lineages and identify the markers."
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed

def step_no_lineage(agent, response_parsed, marker_str):
    """
    Step: If there is no lineage marker found, feed more differential genes.
    """
    query = (
        f"In the updated differential gene list, {marker_str},"
        "can you identify the cell lineages? "
        "If yes, please provide the cell lineages and their markers. "
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed

def step_decide_one_lineage(agent, response_parsed, marker_str):
    """
    Step: If having cell lineages, decide the most likely cell lineage and its markers.
    Args:
        agent: The AI agent for annotation.
        response_parsed: Parsed response from the agent.
        marker_str (str): Comma-separated string of marker genes.
    """
    lineages = ", ".join(response_parsed.cell_lineage) # convert list of text
    query = (
        f"Given a list of potential cell lineage {lineages} and the differential genes {marker_str}, "
        "can you identify the most likely cell lineage? "
        "If yes, please provide the one lineage name and its markers. "
        "If no, ask for more differential genes."
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed



def step_celltypes(agent, response_parsed, marker_str):
    """
    Step: Identify cell types based on the marker genes and the lineage.
    Args:
        agent: The AI agent for annotation.
        response_parsed: Parsed response from the agent.
        marker_str (str): Comma-separated string of marker genes.
    """
    lineage = response_parsed.cell_lineage[0] # get the first lineage
    query = (
        f"Given previously identified cell lineage {lineage}, "
        f"and an enriched differential gene list {marker_str}, can you decide the cell types? "
        "If yes, please provide the cell type names and their markers. "
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed

def step_decide_one_celltype(agent, response_parsed, marker_str):
    """
    Step: If having cell types, decide the most likely cell type and its markers.
    Args:
        agent: The AI agent for annotation.
        response_parsed: Parsed response from the agent.
        marker_str (str): Comma-separated string of marker genes.
    """
    celltypes = ", ".join(response_parsed.cell_type) # convert list of text
    query = (
        f"Given a list of potential cell types {celltypes} and the differential genes {marker_str}, "
        "can you identify the most likely cell type? "
        "If yes, please provide the one cell type name and its markers. "
        "If no, ask for more differential genes."
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed

def step_fcn_anno(agent, response_parsed, marker_str):
    """
    Step: provide functional annotations to the cell type.
    Args:
        agent: The AI agent for annotation.
        response_parsed: Parsed response from the agent.
        marker_str (str): Comma-separated string of marker genes.
    """
    celltype = response_parsed.cell_type[0] # get the first cell type
    query = (
        f"Given the previously identified {celltype} "
        f"and more differential gene list {marker_str}, "
        "can you provide more detailed annotation based on cell development stage, specific functions, activation, "
        "cell cycling, potential response to stimuli, or disease states to this cell type? "
        "Please relate the functions with the genes. "
        "If yes, please provide the refined cell type and the corresponding markers."
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed

def step_decide_one_fcn_anno(agent, response_parsed, marker_str):
    """
    Step: If having functional annotations, decide the most likely functional annotation and its markers.
    Args:
        agent: The AI agent for annotation.
        response_parsed: Parsed response from the agent.
        marker_str (str): Comma-separated string of marker genes.
    """
    celltypes = ", ".join(response_parsed.cell_type) # convert list of text
    query = (
        f"Given the previously identified {celltypes} "
        f"and the differential gene list {marker_str}, "
        "can you identify the most likely cell type? "
        "If yes, please provide the one cell type name and its markers. "
        "If no, ask for more differential genes."
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)
    return response_parsed

### ===========================================================================
### steps used for grouping clusters
def step_group_clusters(agent, clusterID_anno_dict_ls):
    """
    Step: Group clusters based on their annotations.
    Args:
        agent: The AI agent for annotation.
        cluster_anno_dict: Dictionary containing the annotation information.
            - keys are cluster indexes
            - values are annotations
        target_idx: Index of the target cluster.
        query_idx: Index of the query cluster.
    """
    input_list = "\n".join(f"- {entry['id']}: {entry['cell_type']}" for entry in clusterID_anno_dict_ls)

    query = f"Here are the entries:\n{input_list}"
    
    response = agent.run(query)
    response_parsed = agent.parse_output(response)

    return response_parsed

### ===========================================================================
### steps used for cross cluster annotations
def step_cross_refine_single(agent, input_dict, target_idx):
    """
    Step: Refine the annotation for the target cluster
    Args:
        agent: The AI agent for annotation.
        input_dict: Dictionary containing the annotation information.
            input_dict[target_idx]["annotation"]: Annotation of the target cluster.
            input_dict[target_idx]["marker_list"]: List of markers for the target cluster.
            input_dict[target_idx]["gene_list"]: List of genes for the target cluster.
        target_idx: Index of the target cluster.
    """
    
    # extract tagrted cluster annotation and gene list
    target_idx_anno = input_dict[target_idx]["annotation"]
    # markers previously identified
    target_idx_markers = input_dict[target_idx]["marker_list"]
    # new differential genes identified within a group of clusters
    target_idx_genes = input_dict[target_idx]["gene_list"]

    query = (
        f"This cluster is annotated as **{target_idx_anno}** based on gene markers {target_idx_markers}.\n\n "
        "By comparing this cell cluster with other similar cell clusters, " 
        f"there are extra differential genes, which are {target_idx_genes}.\n\n "
        "Can you provide more detailed functional annotations to this cell cluster based on "
        "cell development stage, activation, cell cycling, potential response to stimuli? "
        "If yes, please renew the annotation, provide explanations/analysis, and add new markers. "
        "If you do not change the original annotation, respond no and give explanations/analysis. "
    )
    response = agent.run(query)
    response_parsed = agent.parse_output(response)

    return response_parsed

""" 
### unused function for cross-annotation and merging clusters

# If the response is from step_cross_refine_single is no, perform pairwise comparison
def step_cross_refine_pairwise(agent, input_dict, target_idx, query_idx):

    # Step: cross comparison to decide whether to merge two clusters.

    
    # extract tagrted cluster annotation and gene list
    target_idx_anno = input_dict[target_idx]["annotation"]
    # markers previously identified
    target_idx_markers = input_dict[target_idx]["marker_list"]
    # new differential genes identified within a group of clusters
    target_idx_genes = input_dict[target_idx]["gene_list"]

    # extract the query cluster annotation and gene list
    query_idx_anno = input_dict[query_idx]["annotation"]
    # markers previously identified
    query_idx_markers = input_dict[query_idx]["marker_list"]
    # new differential genes identified within a group of clusters  
    query_idx_genes = input_dict[query_idx]["gene_list"]

    query = (
        "There are two annotated cell clusters in the following. "
        "Compare them and answer questions.\n\n"
        f"Target cell type: {target_idx_anno}.\n"
        f"---Target cell type gene markers: {target_idx_markers}.\n"
        f"---Target cell type extra differential genes: {target_idx_genes}.\n\n"
        f"Query cell type: {query_idx_anno}.\n"
        f"---Query cell type gene markers: {query_idx_markers}.\n"
        f"---Query cell type extra differential genes: {query_idx_genes}.\n\n"
        "Question: should we merge these two clusters/cell types? "
        "If Yes, please provide explanations/analysis and a new cell type for the merged with added new markers. "
        "If No, please provide explanations/analysis."
    )
    
    response = agent.run(query)
    response_parsed = agent.parse_output(response)

    return response_parsed
"""