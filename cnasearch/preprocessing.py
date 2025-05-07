"""
Preprocessing utilities for single-cell data before CNA detection.
"""

import scanpy as sc
import anndata as ad
from typing import List, Optional
from scipy.sparse import issparse


def preprocess(
    adata: sc.AnnData,
    
    # QC parameters
    min_genes_by_counts: int = 200,
    min_counts: int = 1000,
    max_counts: Optional[int] = None,
    max_mito_pct: float = 5.0,
    min_cells_per_gene: Optional[int] = 3,
    
    # gene filtering
    allowed_chromosomes: Optional[List[str]] = None,
    
    # normalization
    target_sum: float = 1e4,
    log1p: bool = True,
    densify: bool = False

    # gene annotation
    annotate_if_missing: bool = True,
    gene_col: Optional[str] = None
) -> sc.AnnData:
    """
    Memory-efficient preprocessing pipeline for scRNA-seq before CNA calling.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object with raw counts
    min_genes_by_counts : int, default=200
        Minimum number of genes per cell
    min_counts : int, default=1000
        Minimum number of counts per cell
    max_counts : Optional[int], default=None
        Maximum number of counts per cell, if None, no upper limit
    max_mito_pct : float, default=5.0
        Maximum percentage of mitochondrial counts per cell
    min_cells_per_gene : Optional[int], default=3
        Minimum number of cells expressing each gene
    allowed_chromosomes : Optional[List[str]], default=None
        Restrict to only these chromosomes, if None, use all
    target_sum : float, default=1e4
        Target sum for library size normalization
    log1p : bool, default=True
        Whether to log1p transform after normalization
    densify : bool, default=False
        Whether to convert sparse matrix to dense
    annotate_if_missing : bool, default=True
        Whether to annotate genes if genomic coordinates are missing
    gene_col : Optional[str], default=None
        Column in adata.var to use for gene IDs when annotating

    Returns
    -------
    AnnData
        Preprocessed AnnData object
    
    Notes
    -----
    - Retains sparse matrix format unless `densify=True`.
    - Uses Scanpy built-ins for filtering and normalization.
    - Requires chromosome, start, end annotations in adata.var
    """
    # 1. Compute QC metrics
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
    )
    
    # 2. Filter cells by mitochondrial content, total counts, and gene counts
    mask = (
        (adata.obs['pct_counts_mt'] < max_mito_pct) &
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['n_genes_by_counts'] >= min_genes_by_counts)
    )
    if max_counts is not None:
        mask &= (adata.obs['total_counts'] <= max_counts)
    adata = adata[mask].copy()
    print(f"[QC] Filtered to {adata.n_obs} cells")
    
    # 3. Filter genes by minimum cells
    if min_cells_per_gene is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        print(f"[QC] Retained {adata.n_vars} genes after min_cells_per_gene={min_cells_per_gene}")
    
    # 4. Annotate, drop genes lacking coordinates, and restrict chromosomes
    coord_cols = ['chromosome', 'start', 'end']
    missing_cols = [col for col in coord_cols if col not in adata.var.columns]
    
    if missing_cols:
        print(f"[QC] Missing genomic coordinate columns: {', '.join(missing_cols)}")
        if annotate_if_missing:
            print("[QC] Annotating genes with genomic coordinates...")
            from .genomics import annotate_genes_mygene
            adata = annotate_genes_mygene(adata, gene_col=gene_col, drop_missing=True)
        else:
            print("[QC] Skipping annotation as annotate_if_missing=False")
            # Initialize empty columns to avoid errors in the next steps
            for col in missing_cols:
                adata.var[col] = None
    
    var = adata.var.dropna(subset=['chromosome', 'start', 'end'])
    if allowed_chromosomes is not None:
        var = var[var['chromosome'].astype(str).isin(allowed_chromosomes)]
    adata = adata[:, var.index].copy()
    print(f"[QC] Retained {adata.n_vars} genes with genomic coords")
    
    # 5. Normalize counts (TP10K)
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
    print(f"[Norm] Scaled to {target_sum} counts per cell")
    
    # 6. Log1p transform
    if log1p:
        sc.pp.log1p(adata)
        print("[Norm] Applied log1p")
    
    # 7. Optional densification
    if densify:
        adata.X = adata.X.toarray() if issparse(adata.X) else adata.X
        print("[Post] Converted to dense array")
    
    print(f"[Preprocess] Completed: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata
