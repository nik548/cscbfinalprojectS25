import numpy as np
from typing import Optional

def log_decay_smoothing(
    B: np.ndarray,
    bin_centers: np.ndarray,
    decay_scale: float = 1e6,
    radius: Optional[int] = 10
) -> np.ndarray:
    """
    Smooth bin-level signals by weighting neighboring bins with a log-distance decay.

    Parameters
    ----------
    B
        Raw bin-expression matrix (cells × bins).
    bin_centers
        Genomic midpoint for each bin.
    decay_scale
        Scale parameter controlling decay steepness.
    radius
        Maximum bin-distance (in bin units) to include; beyond this weight=0.

    Returns
    -------
    S
        Smoothed bin-expression matrix (cells × bins).
    """
    n_bins = B.shape[1]
    # Pairwise genomic distances between bin centers
    D = np.abs(bin_centers[:, None] - bin_centers[None, :])
    # Initialize weight matrix
    W = np.zeros_like(D, dtype=float)
    # Compute log-decay weights for non-zero distances
    mask = D > 0
    W[mask] = 1.0 / np.log1p(D[mask] / decay_scale)
    # Self-weight = 1
    np.fill_diagonal(W, 1.0)
    # Enforce radius cutoff if provided
    if radius is not None:
        idx = np.arange(n_bins)
        M = np.abs(idx[:, None] - idx[None, :]) > radius
        W[M] = 0.0
    # Normalize rows to sum to 1
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    W = W / row_sums
    # Apply smoothing: each cell's bin vector times weight matrix
    return B.dot(W.T)

