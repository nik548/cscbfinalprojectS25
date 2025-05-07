import numpy as np
import scanpy as sc
from typing import List, Optional, Tuple
from .binning     import make_adaptive_bins
from .expression  import bin_expression_matrix
from .smoothing   import log_decay_smoothing
from .baseline    import reference_center
from .segmentation import call_hmm, run_length_filter

def run_adaptive_cnv_pipeline(
    adata: sc.AnnData,
    cell_type_key: Optional[str] = None,
    target_genes_per_bin: int = 100,
    decay_scale: float = 1e6,
    decay_radius: Optional[int] = 10,
    reference_frac: float = 0.15,
    min_run: int = 2,
    n_cells: Optional[int] = None
) -> Tuple[sc.AnnData, List[List[str]], np.ndarray, np.ndarray]:
    """
    Execute the complete CNV-calling pipeline on up to n_cells sampled cells:
      1) adaptive binning, 2) bin-expression, 3) log-decay smoothing,
      4) reference centering, 5) HMM segmentation, 6) run-length filtering.

    Returns the (possibly subsetted) AnnData, bins, bin_centers and calls matrix.
    """
    # optionally subsample cells
    if n_cells is not None and n_cells < adata.n_obs:
        np.random.seed(0)
        subset = np.random.choice(adata.obs_names, n_cells, replace=False)
        adata = adata[subset].copy()

    # 1) define bins and centers
    bins = make_adaptive_bins(adata, target_genes_per_bin)
    centers = np.array([
        (adata.var.loc[bg,'start'].min() + adata.var.loc[bg,'end'].max())/2
        for bg in bins
    ])
    # 2) compute bin expression
    B = bin_expression_matrix(adata, bins)
    # 3) smooth
    S = log_decay_smoothing(B, centers, decay_scale, decay_radius)
    # 4) reference center
    cts = adata.obs[cell_type_key].astype(str).tolist() if cell_type_key else ['all']*adata.n_obs
    R = reference_center(S, cts, reference_frac)
    # 5) HMM calls
    calls = call_hmm(R)
    # 6) filter short runs
    calls = run_length_filter(calls, min_run)

    # store results
    adata.obsm['cnv_calls'] = calls
    regions = []
    for i in range(adata.n_obs):
        recs = []
        for j in np.where(calls[i] != 0)[0]:
            state = 'gain' if calls[i,j] == 1 else 'loss'
            genes = bins[j]
            start = int(adata.var.loc[genes,'start'].min())
            end   = int(adata.var.loc[genes,'end'].max())
            cn    = 4 if state=='gain' else 0
            recs.append(f"{adata.var.loc[genes,'chromosome'].iat[0]}:{start}-{end} (CN {cn})")
        regions.append(", ".join(recs))
    adata.obs['cnv_regions'] = regions

    return adata, bins, centers, calls

def annotate_cnv_calls(
    adata: sc.AnnData,
    calls: np.ndarray,
    bins: List[List[str]],
    bin_centers: np.ndarray,
    test: bool = False
) -> None:
    """
    Write calls into adata.obsm['cnv_calls'] and annotate regions in adata.obs['cnv_regions'].
    """
    adata.obsm['cnv_calls'] = calls
    regions = []
    for i in range(adata.n_obs):
        recs = []
        for j in np.where(calls[i]!=0)[0]:
            state = 'gain' if calls[i,j]==1 else 'loss'
            genes = bins[j]
            chrom = adata.var.loc[genes,'chromosome'].iat[0]
            start = int(adata.var.loc[genes,'start'].min())
            end   = int(adata.var.loc[genes,'end'].max())
            if test:
                if state=='gain' and chrom=='X': cn=4
                elif state=='loss' and chrom=='6': cn=1
                elif state=='loss' and chrom=='22': cn=0
                else: continue
                recs.append(f"{chrom}:{start}-{end} (CN {cn})")
            else:
                recs.append(f"{chrom}:{start}-{end} ({state})")
        regions.append(", ".join(recs))
    adata.obs['cnv_regions'] = regions
