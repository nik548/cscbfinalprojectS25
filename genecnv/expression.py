import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from typing import List

def bin_expression_matrix(
    adata: sc.AnnData,
    bins: List[List[str]]
) -> np.ndarray:
    """
    Compute the mean expression per cell within each gene bin.

    Parameters
    ----------
    adata
        Annotated data matrix with expression in .X.
    bins
        List of gene bins (from make_adaptive_bins).

    Returns
    -------
    B
        Array of shape (n_cells, n_bins) of average expression.
    """
    # Extract dense matrix if sparse
    X = adata.X.toarray() if issparse(adata.X) else adata.X
    # Map gene name to column index
    idx = {g:i for i,g in enumerate(adata.var_names)}
    # Initialize output
    B = np.zeros((adata.n_obs, len(bins)), float)
    # Compute mean across genes in each bin
    for j, bin_genes in enumerate(bins):
        cols = [idx[g] for g in bin_genes]
        B[:, j] = X[:, cols].mean(axis=1)
    return B
