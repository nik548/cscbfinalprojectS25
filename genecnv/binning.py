import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from typing import List, Optional, Dict, Tuple
from hmmlearn import hmm

# ──────────────────────────────────────────────────────────────────────────────
# 1) Adaptive gene‑centric bins
# ──────────────────────────────────────────────────────────────────────────────
def make_adaptive_bins(
    adata: sc.AnnData,
    target_genes_per_bin: int = 100
) -> List[List[str]]:
    """
    Partition genes into sequential bins of approximately equal gene count,
    ordering by chromosome and genomic start position.

    Parameters
    ----------
    adata
        Annotated data matrix with .var['chromosome'] and .var['start'] fields.
    target_genes_per_bin
        Approximate number of genes per adaptive bin.

    Returns
    -------
    bins
        List of bins, each a list of gene names.
    """
    # Filter to genes with valid chromosome/start information
    df = adata.var.dropna(subset=['chromosome','start']).copy()
    # Create numeric ordering for chromosomes (1-22, X=23, Y=24, others=25)
    df['chrom_ord'] = df['chromosome'].map(
        lambda c: int(c) if str(c).isdigit() else {'X':23,'Y':24}.get(c,25)
    )
    # Sort by chromosome order then by genomic coordinate
    df = df.sort_values(['chrom_ord','start'])
    genes = df.index.tolist()
    # Chunk gene list into bins of specified size
    return [genes[i:i+target_genes_per_bin] for i in range(0, len(genes), target_genes_per_bin)]
