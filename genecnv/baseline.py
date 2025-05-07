import numpy as np
from typing import List

def reference_center(
    S: np.ndarray,
    cell_types: List[str],
    reference_frac: float = 0.15
) -> np.ndarray:
    """
    Subtract per-cell-type reference mean (diploid baseline) from smoothed signals.

    Parameters
    ----------
    S
        Smoothed bin-expression (cells × bins).
    cell_types
        List of cell-type labels for each cell.
    reference_frac
        Fraction of lowest-variance cells per type to use as reference.

    Returns
    -------
    R
        Centered signal (cells × bins).
    """
    R = np.zeros_like(S)
    # Process each cell type separately
    for ct in np.unique(cell_types):
        ids = np.where(np.array(cell_types) == ct)[0]
        # Compute variance across bins for each cell
        var = S[ids].var(axis=1)
        # Select lowest-variance cells as reference
        n_ref = max(int(len(ids) * reference_frac), 1)
        ref = ids[np.argsort(var)[:n_ref]]
        # Compute mean reference profile
        mu = S[ref].mean(axis=0)
        # Subtract baseline
        R[ids] = S[ids] - mu
    return R
