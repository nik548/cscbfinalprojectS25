"""
Hidden Markov Model (HMM) utilities for CNA detection.
"""

import numpy as np
from typing import List, Dict, Tuple
from hmmlearn import hmm


def smooth_calls(calls: np.ndarray, min_run: int = 2) -> np.ndarray:
    """
    Filter CNA calls by requiring a minimum run length.
    
    Parameters
    ----------
    calls : ndarray
        Array of shape (n_cells, n_windows) with values {-1, 0, +1}
    min_run : int, default=2
        Minimum number of consecutive windows required to keep a call
    
    Returns
    -------
    out : ndarray
        Filtered calls with short runs removed
    """
    out = np.zeros_like(calls)
    for i in range(calls.shape[0]):
        for state in (+1, -1):
            idx = np.where(calls[i] == state)[0]
            if idx.size == 0: 
                continue
            runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
            for run in runs:
                if len(run) >= min_run:
                    out[i, run] = state
    return out


def call_cnas_hmm_viterbi_by_chrom(
    window_expr: np.ndarray,
    reference_idx: np.ndarray,
    window_labels: List[str],
    chrom_models: Dict[str, hmm.GaussianHMM] = None,
    n_components: int = 3,
    cov_type: str = 'diag'
) -> np.ndarray:
    """
    Call CNAs using HMM with Viterbi decoding, fitting one model per chromosome.
    
    Parameters
    ----------
    window_expr : ndarray
        Array of shape (n_cells, n_windows) with window expression values
    reference_idx : ndarray
        Indices of reference cells
    window_labels : List[str]
        List of window identifiers (chr_window)
    chrom_models : Dict[str, hmm.GaussianHMM], default=None
        Pre-fitted HMM models per chromosome, if None will fit new models
    n_components : int, default=3
        Number of HMM states (usually 3: loss, neutral, gain)
    cov_type : str, default='diag'
        Covariance type for HMM
    
    Returns
    -------
    calls : ndarray
        Array of shape (n_cells, n_windows) with calls {-1, 0, +1}
    
    Notes
    -----
    - Fits one HMM per chromosome
    - Assigns states to loss/neutral/gain based on ordering of means
    """
    # Assemble Z-scored expression matrix
    Z = window_expr.astype(np.float32).copy()
    
    # Z-score globally by reference
    mu = Z[reference_idx].mean(axis=0)
    sigma = Z[reference_idx].std(axis=0, ddof=1)
    sigma[sigma == 0] = 1.0
    Z = (Z - mu) / sigma
    
    # Group window indices by chromosome
    chrom_to_idxs = {}
    for j, lab in enumerate(window_labels):
        chrom = lab.split('_', 1)[0]
        chrom_to_idxs.setdefault(chrom, []).append(j)
    
    calls = np.zeros_like(Z, dtype=int)
    
    # Fit & decode per chromosome
    for chrom, idxs in chrom_to_idxs.items():
        # Use provided model or train new HMM
        model = hmm.GaussianHMM(
            n_components=n_components,
            covariance_type=cov_type,
            n_iter=20, 
            random_state=0
        )
        
        # Reshape data for training (n_cells*n_windows_chrom, 1)
        data = Z[:, idxs].reshape(-1, 1)
        model.fit(data)
        
        # Decode each cell
        for i in range(Z.shape[0]):
            path = model.predict(Z[i, idxs].reshape(-1, 1))
            
            # Assign states based on ordering of means
            order = np.argsort(model.means_.flatten())
            loss_s, gain_s = order[0], order[-1]
            
            # Map to calls
            calls[i, np.array(idxs)[path == gain_s]] = +1
            calls[i, np.array(idxs)[path == loss_s]] = -1
    
    return calls


def fit_and_get_posteriors(
    W: np.ndarray,
    labels: List[str],
    reference_idx: np.ndarray,
    n_components: int = 3,
    cov_type: str = 'diag'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit HMMs and get posterior probabilities for gain and loss states.
    
    Parameters
    ----------
    W : ndarray
        Array of shape (n_cells, n_windows) with window expression values
    labels : List[str]
        List of window identifiers (chr_window)
    reference_idx : ndarray
        Indices of reference cells
    n_components : int, default=3
        Number of HMM states (usually 3: loss, neutral, gain)
    cov_type : str, default='diag'
        Covariance type for HMM
    
    Returns
    -------
    post_gain : ndarray
        Array of shape (n_cells, n_windows) with gain state posteriors
    post_loss : ndarray
        Array of shape (n_cells, n_windows) with loss state posteriors
    """
    # Z-score by reference
    Z = W.astype(np.float32).copy()
    mu = Z[reference_idx].mean(axis=0)
    sigma = Z[reference_idx].std(axis=0, ddof=1)
    sigma[sigma == 0] = 1.0
    Z = (Z - mu) / sigma
    
    # Group windows by chromosome
    chrom_to_idxs: Dict[str, List[int]] = {}
    for j, lab in enumerate(labels):
        chrom = lab.split('_', 1)[0]
        chrom_to_idxs.setdefault(chrom, []).append(j)
    
    post_gain = np.zeros_like(Z)
    post_loss = np.zeros_like(Z)
    
    # Fit per chromosome
    for chrom, idxs in chrom_to_idxs.items():
        data = Z[:, idxs].reshape(-1, 1)
        model = hmm.GaussianHMM(
            n_components=n_components,
            covariance_type=cov_type,
            n_iter=20, 
            random_state=0
        )
        model.fit(data)
        
        # Identify loss/gain states by means
        order = np.argsort(model.means_.flatten())
        loss_s, gain_s = order[0], order[-1]
        
        # Compute posterior on full data
        post = model.predict_proba(data).reshape(Z.shape[0], len(idxs), n_components)
        post_gain[:, idxs] = post[:, :, gain_s]
        post_loss[:, idxs] = post[:, :, loss_s]
    
    return post_gain, post_loss