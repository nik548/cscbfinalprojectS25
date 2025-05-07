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
) -> sc.AnnData:
    """
    Memory-efficient preprocessing pipeline for scRNA-seq before CNA calling.

    - Retains sparse matrix format unless `densify=True`.
    - Uses Scanpy built-ins for filtering and normalization.
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

    # 4. Drop genes lacking coordinates and restrict chromosomes
    var = adata.var.dropna(subset=['chromosome','start','end'])
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
