import numpy as np
from hmmlearn import hmm

def call_hmm(R: np.ndarray) -> np.ndarray:
    """
    Fit a 3-state Gaussian HMM across all cells×bins and decode gain/loss states per cell.

    Parameters
    ----------
    R
        Reference-centered signals (cells × bins).

    Returns
    -------
    calls
        Integer array (cells × bins) with -1=loss, 0=neutral, +1=gain.
    """
    # Replace NaNs with zero
    R = np.nan_to_num(R, nan=0.0)
    # Fit Gaussian HMM on flattened data
    model = hmm.GaussianHMM(n_components=3, covariance_type='diag', n_iter=20)
    model.fit(R.reshape(-1, 1))
    # Identify which state corresponds to loss vs gain by ordering means
    order = np.argsort(model.means_.flatten())
    loss_s, gain_s = order[0], order[-1]
    # Initialize call matrix
    calls = np.zeros_like(R, int)
    # Decode each cell's Viterbi path
    for i in range(R.shape[0]):
        path = model.predict(R[i].reshape(-1,1))
        calls[i, path == gain_s] = +1
        calls[i, path == loss_s] = -1
    return calls

def run_length_filter(calls: np.ndarray, min_run: int = 2) -> np.ndarray:
    """
    Remove isolated CNA calls shorter than min_run consecutive bins.

    Parameters
    ----------
    calls
        Raw calls (cells × bins) from call_hmm.
    min_run
        Minimum consecutive bins to retain a call region.

    Returns
    -------
    filtered
        Filtered calls array (cells × bins).
    """
    out = np.zeros_like(calls)
    for i in range(calls.shape[0]):
        for state in (+1, -1):
            idxs = np.where(calls[i] == state)[0]
            runs = np.split(idxs, np.where(np.diff(idxs) > 1)[0] + 1)
            for run in runs:
                if len(run) >= min_run:
                    out[i, run] = state
    return out
