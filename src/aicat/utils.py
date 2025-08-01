import scanpy as sc
import pandas as pd
import subprocess

def normalize_count_data(data, n_top_genes=2000):
    """
    Calculate differential genes from the given data file.

    Args:
        data (AnnData): AnnData object containing the single-cell RNA-seq data.
        n_top_genes (int): Number of top genes to keep.

    Returns:
        dict: Dictionary with cluster names as keys and lists of top genes as values.
    """
    # Preprocess the data
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=n_top_genes)
    
    # Focus on the highly variable genes
    data_hvg = data[:, data.var["highly_variable"]].copy()
    data_hvg.raw = None

    return data_hvg

def sub_clustering(data_hvg, cell_type_col, cluster_idx,
                   key_added="subcluster",
                   resolution=0.8):
    """
    Perform sub-clustering on the given data.

    Args:
        data_hvg (AnnData): Preprocessed AnnData object.
        cell_type_col (str): Column name representing cell types.
        cluster_idx (int): Index of the cluster to be sub-clustered.
        resolution (float): Resolution parameter for clustering. Higher, more clusters.

    Returns:
        AnnData: Sub-clustered AnnData object.
    """
    # Subset the data for the specified cluster
    data_sub = data_hvg[data_hvg.obs[cell_type_col] == cluster_idx].copy()
    
    # Perform PCA and UMAP
    sc.pp.pca(data_sub)
    sc.pp.neighbors(data_sub)
    sc.tl.umap(data_sub) # Do UMAP incase of visualization
    sc.tl.leiden(data_sub, key_added=key_added, resolution=resolution,
                 flavor="igraph",
                 directed=False)  # Use igraph for Leiden clustering

    return data_sub

def calculate_diff_genes(data_hvg, cell_type_col,
                         pval_cutoff=0.05, logfc_cutoff=1,
                         n_genes_to_keep=100,
                         save_path=None):
    """
    Calculate differential genes for the given data.

    Args:
        data_hvg (AnnData): Preprocessed AnnData object.
        cell_type_col (str): Column name representing cell types.
        pval_cutoff (float): Adjusted p-value cutoff.
        logfc_cutoff (float): Log fold change cutoff.
        n_genes_to_keep (int): Number of top genes to keep for each cluster, sorted by logfoldchange.

    Returns:
        dict: Dictionary with cluster names as keys and lists of top genes as values.
    """
    # Obtain cluster-specific differentially expressed genes
    sc.tl.rank_genes_groups(data_hvg, groupby=cell_type_col, method="wilcoxon")
    
    # Extract results for all clusters
    df = sc.get.rank_genes_groups_df(data_hvg, None)
    
    # Filter for adjusted p-value < 0.05 and log fold change > 1
    filtered_df = df[(df["pvals_adj"] < pval_cutoff) & (df["logfoldchanges"] > logfc_cutoff)]

    # Create a dictionary with cluster names as keys and lists of top genes as values
    top_genes_dict = (
        filtered_df.sort_values(by="logfoldchanges", ascending=False)
          .groupby("group", observed=True)
          .head(n_genes_to_keep)
          .groupby("group", observed=True)["names"]
          .apply(list)
          .to_dict()
    )

    # Save the results to a JSON file if a save path is provided
    if save_path:
        import json
        with open(save_path, "w") as f:
            json.dump(top_genes_dict, f, indent=2)


    return top_genes_dict

def calculate_within_group_markers(data_hvg, cell_type_col,
                                   grouped_cluster_idx_dict,
                                   pval_cutoff=0.05, logfc_cutoff=1,
                                   n_genes_to_keep=100,
                                   save_path=None):
    """
    Calculate differential genes for each group of clusters.

    Args:
        data_hvg (AnnData): Preprocessed AnnData object.
        cell_type_col (str): Column name representing cell types.
        grouped_cluster_idx_dict (dict): Dictionary with cluster indices as keys, groups as values.
        pval_cutoff (float): Adjusted p-value cutoff.
        logfc_cutoff (float): Log fold change cutoff.
        n_genes_to_keep (int): Number of top genes to keep for each group, sorted by logfoldchange.
        save_path (str): Path to save the results.
    
    Returns:
        dict: Dictionary with group names as keys and lists of top genes as values.
    """
    
    within_group_markers_dict = {}
    for group, cluster_idx_list in grouped_cluster_idx_dict.items():
        if len(cluster_idx_list) < 2:
            continue # Skip groups with fewer than 2 clusters

        # Subset the data for the specified clusters
        data_sub = data_hvg[data_hvg.obs[cell_type_col].isin(cluster_idx_list)].copy()
        # Keep only highly variable genes used in the original adata
        data_sub = data_sub[:, data_sub.var.highly_variable].copy()

        within_group_markers = calculate_diff_genes(data_sub, cell_type_col,
                                                    pval_cutoff=pval_cutoff, 
                                                    logfc_cutoff=logfc_cutoff,
                                                    n_genes_to_keep=n_genes_to_keep,
                                                    save_path=None)
        within_group_markers_dict[group] = within_group_markers
    
    # Save the results to a JSON file if a save path is provided
    if save_path:
        import json
        with open(save_path, "w") as f:
            json.dump(within_group_markers_dict, f, indent=2)

    return within_group_markers_dict



def chat_to_markdown(chat_history):
    """
    Convert chat history to markdown format.
    """
    md_lines = []
    for role, message in chat_history:
        role_title = "**USER:**" if role.lower() == "user" else "**AI:**"
        md_lines.append(f"{role_title}\n{message.strip()}\n")
    return "\n".join(md_lines)

def convert_md_to_pdf_with_pandoc(md_file, pdf_file):
    """
    Convert markdown file to PDF using Pandoc.

    **Requirements: pandoc**
    •	macOS: brew install pandoc
	•	Ubuntu: sudo apt install pandoc
	•	Windows: Download from pandoc.org
    """
    subprocess.run(["pandoc", md_file, "-o", pdf_file], check=True)

def collect_annotation_single(single_res):
    """
    Collect annotation results from a single cell type dictionary.

    Parameters:
    -----------
    single_res : dict
        Dictionary where each key is a ground truth label (e.g., cluster or cell type),
        and each value is a nested dictionary with annotation categories:
            - 'lineage' -> {'cell_lineage': str or [str]}
            - 'celltype' -> {'cell_type': str or [str]}
            - 'fcn_anno' -> {'cell_type': str or [str]}

    Returns:
    --------
    pd.DataFrame
        A DataFrame with columns: 'Ground_truth', 'lineage', 'celltype', 'fcn_anno'.
    """
    def extract_value(val):
        # If it's a list, return the first element; else return the value directly
        if isinstance(val, list):
            return val[0] if val else "NA"
        return val if val is not None else "NA"

    results = []

    for ground_truth, annotations in single_res.items():
        lineage = extract_value(annotations.get("lineage", {}).get("cell_lineage"))
        celltype = extract_value(annotations.get("celltype", {}).get("cell_type"))
        fcn = extract_value(annotations.get("fcn_anno", {}).get("cell_type"))
        summary = extract_value(annotations.get("fcn_anno", {}).get("summary"))

        results.append([ground_truth, lineage, celltype, fcn, summary])

    return pd.DataFrame(results, columns=["Ground_truth", "lineage", 
                                          "celltype", "fcn_anno", "summary"])

def format_query_prompt(row) -> str:
    """
    Formats a natural language query prompt for the annotation evaluation agent
    based on available fields in a row.
    
    Parameters:
    - row: A pandas Series with the following columns (some may be missing or None):
        - 'ground_truth'
        - 'lineage'
        - 'cell_type'
        - 'functional'
        - 'summary'

    Returns:
    - A formatted string prompt
    """
    ground_truth = row.get('Ground_truth', 'None')
    lineage = row.get('lineage', 'None') or 'None'
    cell_type = row.get('celltype', 'None') or 'None'
    functional = row.get('fcn_anno', 'None') or 'None'
    summary = row.get('summary', 'None') or 'None'

    prompt = f"""
        Ground truth: {ground_truth}
        Lineage annotation: {lineage}
        Cell type annotation: {cell_type}
        Functional annotation: {functional}
        Annotation summary: {summary}

        Please evaluate the above annotations using the defined scoring system. Return:
        - lineage_score (if lineage is provided)
        - cell_type_score (if cell type is provided)
        - functional_score (if functional annotation is provided)
        - final_score
        - justification for each score
        - sources (if any)
        - tools_used (if any)
        """.strip()

    return prompt