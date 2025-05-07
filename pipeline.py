"""
Full pipeline implementation for CNAsearch.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import List, Dict, Union, Optional

from .windows import window_genes_by_distance, compute_window_expression, select_reference_cells
from .hmm import fit_and_get_posteriors, smooth_calls


def run_cna_pipeline_robust(
    adata: sc.AnnData,
    window_size: int = 1_000_000,
    min_genes_per_window: int = 10,
    ref_method: str = 'kmeans',
    ref_frac: float = 0.15,
    min_run: int = 2,
    prob_thresh: float = 0.5
) -> List[Dict]:
    """
    Run the full CNAsearch pipeline with probabilistic HMM-based calling.
    
    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object
    window_size : int, default=1_000_000
        Size of each genomic window in base pairs
    min_genes_per_window : int, default=10
        Minimum number of genes required in a window
    ref_method : str, default='kmeans'
        Method for selecting reference cells, either 'kmeans' or 'variance'
    ref_frac : float, default=0.15
        Fraction of cells to select when using 'variance' method
    min_run : int, default=2
        Minimum number of consecutive windows required to keep a call
    prob_thresh : float, default=0.5
        Probability threshold for calling CNAs from HMM posteriors
    
    Returns
    -------
    summary : List[Dict]
        List of dictionaries with CNA events
    
    Notes
    -----
    Each entry in the summary contains:
    - cell_id: identifier of the cell
    - chromosome: chromosome name
    - window_start: start position of the window
    - window_end: end position of the window
    - call: 'gain' or 'loss'
    """
    # 1) Tile & compute window expression
    var_win, window_bounds = window_genes_by_distance(
        adata, 
        window_size=window_size, 
        min_genes_per_window=min_genes_per_window
    )
    W, labels = compute_window_expression(adata, var_win)
    
    # 2) Pick reference cells
    ref_idx = select_reference_cells(W, method=ref_method, frac=ref_frac)
    
    # 3) Fit HMMs once & get posteriors
    post_gain, post_loss = fit_and_get_posteriors(W, labels, ref_idx)
    
    # 4) Threshold posteriors into calls
    calls = np.zeros_like(post_gain, dtype=int)
    calls[post_gain >= prob_thresh] = +1
    calls[post_loss >= prob_thresh] = -1
    
    # 5) Smooth short runs
    calls = smooth_calls(calls, min_run=min_run)
    
    # 6) Summarize into list of dicts
    summary = []
    for i, cell in enumerate(adata.obs_names):
        for j, lbl in enumerate(labels):
            if calls[i, j] == 0:
                continue
            wb = window_bounds.loc[lbl]
            summary.append({
                'cell_id': cell,
                'chromosome': wb['chromosome'],
                'window_start': int(wb['window_start']),
                'window_end': int(wb['window_end']),
                'call': 'gain' if calls[i, j] == 1 else 'loss'
            })
    
    return summary


def add_cna_annotations_to_obs(
    adata: sc.AnnData,
    cna_df: pd.DataFrame,
    window_size: int = 1_000_000,
    min_genes_per_window: int = 10
) -> sc.AnnData:
    """
    Add CNA annotations to the AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to annotate
    cna_df : DataFrame
        DataFrame with CNA events from run_cna_pipeline_robust
    window_size : int, default=1_000_000
        Size of each genomic window in base pairs
    min_genes_per_window : int, default=10
        Minimum number of genes required in a window
    
    Returns
    -------
    adata : AnnData
        AnnData with CNA annotations in obs
    
    Notes
    -----
    Adds the following columns to adata.obs:
    - cna_regions: comma-separated list of CNA events
    - cna_count: number of CNA events
    - has_cna: boolean indicating presence of any CNA
    """
    # Create a dict of cell_id -> list of CNA regions
    cell_to_cnas = {}
    for _, row in cna_df.iterrows():
        cell = row['cell_id']
        chrom = row['chromosome']
        start = row['window_start']
        end = row['window_end']
        call = row['call']
        
        # Format: chr1:1000000-2000000(gain)
        region = f"{chrom}:{start}-{end}({call})"
        
        cell_to_cnas.setdefault(cell, []).append(region)
    
    # Add to obs
    cna_regions = pd.Series(
        {cell: ",".join(regions) for cell, regions in cell_to_cnas.items()},
        index=adata.obs_names
    ).fillna("")
    
    adata.obs['cna_regions'] = cna_regions
    adata.obs['cna_count'] = cna_regions.str.count(',') + (cna_regions != "").astype(int)
    adata.obs['has_cna'] = adata.obs['cna_count'] > 0
    
    print(f"Added CNA annotations to {adata.obs['has_cna'].sum()} cells")
    return adata
