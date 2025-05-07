"""
Genomic windowing and expression calculation utilities.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse
from typing import List, Dict, Tuple
from sklearn.cluster import KMeans


def window_genes_by_distance(
    adata: sc.AnnData,
    window_size: int = 1_000_000,
    min_genes_per_window: int = 10
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Partition genes into genomic windows based on their coordinates.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with chromosome, start, end in var
    window_size : int, default=1_000_000
        Size of each genomic window in base pairs
    min_genes_per_window : int, default=10
        Minimum number of genes required in a window
    
    Returns
    -------
    var_windowed : DataFrame
        Gene annotations with window assignments
    window_bounds : DataFrame
        DataFrame with window boundary information
    
    Notes
    -----
    This function requires chromosome, start, and end annotations in adata.var
    """
    var = adata.var[['chromosome', 'start', 'end']].dropna().copy()
    var['window_id'] = (var['start'] // window_size).astype(int)
    var['chr_window'] = var['chromosome'].astype(str) + '_' + var['window_id'].astype(str)
    
    # Keep only windows with sufficient genes
    keep = var['chr_window'].value_counts().loc[lambda x: x >= min_genes_per_window].index
    var_windowed = var[var['chr_window'].isin(keep)].copy()
    
    # Calculate window boundaries
    wb = (
        var_windowed.groupby('chr_window')
        .agg(chromosome=('chromosome', 'first'),
             window_id=('window_id', 'first'),
             start_min=('start', 'min'),
             end_max=('end', 'max'))
    )
    wb['window_start'] = wb['window_id'] * window_size
    wb['window_end'] = wb['window_start'] + window_size - 1
    
    # Sort windows by chromosome and position
    def chr_to_ord(c):
        try: 
            return int(c)
        except: 
            return {'X': 23, 'Y': 24}.get(c, 25)
    
    tmp = wb.reset_index()
    tmp[['chrom', 'win']] = tmp['chr_window'].str.split('_', expand=True)
    tmp['chrom_ord'] = tmp['chrom'].map(chr_to_ord)
    tmp['win'] = tmp['win'].astype(int)
    tmp = tmp.sort_values(['chrom_ord', 'win']).drop(['chrom', 'win', 'chrom_ord'], axis=1)
    window_bounds = tmp.set_index('chr_window')
    
    return var_windowed, window_bounds


def compute_window_expression(
    adata: sc.AnnData,
    var_windowed: pd.DataFrame
) -> Tuple[np.ndarray, List[str]]:
    """
    Compute mean expression per genomic window.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with expression data
    var_windowed : DataFrame
        Gene annotations with window assignments from window_genes_by_distance
    
    Returns
    -------
    window_expr : ndarray
        Array of shape (n_cells, n_windows) with window expression values
    window_labels : List[str]
        List of window identifiers
    
    Notes
    -----
    - Uses sparse matrix operations for efficiency
    - Expression of each gene is weighted by 1/n where n is the number of genes in the window
    """
    window_labels = var_windowed['chr_window'].unique().tolist()
    win_idx = {w: i for i, w in enumerate(window_labels)}
    gene_idx = {g: i for i, g in enumerate(adata.var_names)}
    
    # Create a genes × windows matrix for aggregation
    rows, cols, data = [], [], []
    counts = var_windowed['chr_window'].value_counts()
    
    for gene, row in var_windowed.iterrows():
        rows.append(gene_idx[gene])
        cols.append(win_idx[row['chr_window']])
        # Weight by 1/n where n is the number of genes in the window
        data.append(1.0 / counts[row['chr_window']])
    
    G = csr_matrix((data, (rows, cols)), shape=(adata.n_vars, len(window_labels)))
    
    # Matrix multiplication to get window expression
    X = adata.X if issparse(adata.X) else csr_matrix(adata.X)
    W = X.dot(G)
    
    return W.toarray(), window_labels


def select_reference_cells(
    window_expr: np.ndarray,
    method: str = 'kmeans',
    frac: float = 0.15,
    n_clusters: int = 2
) -> np.ndarray:
    """
    Select reference cells for baseline expression.
    
    Parameters
    ----------
    window_expr : ndarray
        Array of shape (n_cells, n_windows) with window expression values
    method : str, default='kmeans'
        Method for selecting reference cells, either 'kmeans' or 'variance'
    frac : float, default=0.15
        Fraction of cells to select when using 'variance' method
    n_clusters : int, default=2
        Number of clusters for kmeans method
    
    Returns
    -------
    ref_idx : ndarray
        Indices of selected reference cells
    
    Notes
    -----
    - 'variance' method selects cells with lowest overall variance
    - 'kmeans' method clusters cells and selects the cluster with lowest variance
    """
    if method == 'variance':
        v = np.nanvar(window_expr, axis=1)
        n = max(int(len(v) * frac), 1)
        return np.argsort(v)[:n]
    
    # KMeans clustering
    km = KMeans(n_clusters=n_clusters, random_state=0).fit(window_expr)
    labs = km.labels_
    
    # Select cluster with lowest variance
    var_by_cluster = [np.nanvar(window_expr[labs == k], axis=1).mean() for k in range(n_clusters)]
    ref_cl = int(np.argmin(var_by_cluster))
    
    return np.where(labs == ref_cl)[0]
