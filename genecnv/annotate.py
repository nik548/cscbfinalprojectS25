"""
Genomic annotation utilities for genecnv.
"""

import anndata as ad
import pandas as pd
from typing import List


def annotate_genes_mygene(
    adata,
    gene_col=None,
    scopes=['ensembl.gene', 'symbol', 'entrezgene'],
    species='human'
):
    """
    Annotate genes with genomic coordinates using MyGene.info API,
    automatically adding start/end columns and filtering to standard chromosomes.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with genes to annotate
    gene_col : str, default=None
        Column in adata.var to use for gene IDs. If None, uses the var index.
    scopes : list, default=['ensembl.gene', 'symbol', 'entrezgene']
        Fields to search for gene IDs (see MyGene.info documentation)
    species : str, default='human'
        Species for MyGene query
    
    Returns
    -------
    AnnData
        Filtered AnnData with chromosome, start, end annotations in var,
        containing only genes with valid annotations on standard chromosomes.
    """
    try:
        import mygene
        import pandas as pd
        import re
    except ImportError:
        raise ImportError("Please install required packages: pip install mygene pandas")
    
    # Create MyGene client
    mg = mygene.MyGeneInfo()
    
    # Create a copy of the input AnnData
    adata_copy = adata.copy()
    
    # Determine which identifiers to use for querying
    if gene_col is None:
        # Use index as gene identifiers
        gene_identifiers = list(adata_copy.var.index)
    else:
        # Use specified column
        if gene_col not in adata_copy.var.columns:
            raise ValueError(f"Column '{gene_col}' not found in adata.var")
        gene_identifiers = adata_copy.var[gene_col].tolist()
    
    # Clean gene identifiers (remove version suffixes and duplicates)
    clean_ids = []
    id_map = {}  # Maps clean IDs to original row indices
    
    print(f"Preparing to query {len(gene_identifiers)} gene identifiers...")
    
    for i, gene_id in enumerate(gene_identifiers):
        if gene_id is None or pd.isna(gene_id):
            continue
            
        # Convert to string and remove version suffixes (e.g., ENSG00000123456.1 -> ENSG00000123456)
        clean_id = re.sub(r'\.\d+$', '', str(gene_id))
        
        # Store mapping from clean ID to row index
        if clean_id not in id_map:
            clean_ids.append(clean_id)
            id_map[clean_id] = []
        id_map[clean_id].append(i)
    
    # Query MyGene.info
    print(f"Querying MyGene.info for {len(clean_ids)} unique gene identifiers...")
    
    results = mg.querymany(
        clean_ids,
        scopes=scopes,
        fields=['chromosome', 'genomic_pos', 'genomic_pos_hg19'],
        species=species,
        returnall=True
    )['out']
    
    # Process results and create annotation DataFrames
    annotations = {}
    
    for result in results:
        # Skip not found
        if 'notfound' in result and result['notfound']:
            continue
            
        query = result.get('query', '')
        
        # Try both genomic_pos fields
        gp = result.get('genomic_pos')
        if gp is None:
            gp = result.get('genomic_pos_hg19')
            
        # Skip if no position info
        if gp is None:
            continue
            
        # Handle different structures of genomic_pos
        if isinstance(gp, list):
            if not gp:
                continue
            gp = gp[0]  # Take first entry
            
        # Skip if not a dictionary
        if not isinstance(gp, dict):
            continue
            
        # Extract coordinates
        chrom = result.get('chromosome') or gp.get('chr')
        start = gp.get('start')
        end = gp.get('end')
        
        # Skip if missing any coordinate
        if None in (chrom, start, end):
            continue
            
        # Store annotations for all occurrences of this gene ID
        if query in id_map:
            for idx in id_map[query]:
                annotations[idx] = {
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'genomic_pos': f"{start}-{end}"
                }
    
    # Define standard chromosomes
    std_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    # Create a new var DataFrame with annotations
    new_var = adata_copy.var.copy()
    
    # Initialize coordinate columns
    for col in ['chromosome', 'start', 'end', 'genomic_pos']:
        new_var[col] = None
    
    # Add annotations
    for idx, annot in annotations.items():
        for col, val in annot.items():
            new_var.iloc[idx, new_var.columns.get_loc(col)] = val
    
    # Filter to standard chromosomes and drop missing annotations
    valid_mask = (
        new_var['chromosome'].isin(std_chroms) & 
        ~new_var['start'].isna() & 
        ~new_var['end'].isna()
    )
    
    # Report statistics
    total_genes = len(new_var)
    annotated_genes = (~new_var['chromosome'].isna()).sum()
    std_chrom_genes = valid_mask.sum()
    
    print(f"Annotation summary:")
    print(f"  Total genes: {total_genes}")
    print(f"  Successfully annotated: {annotated_genes} ({annotated_genes/total_genes:.1%})")
    print(f"  Standard chromosomes (1-22,X,Y): {std_chrom_genes} ({std_chrom_genes/total_genes:.1%})")
    
    # Filter the AnnData object
    filtered_var = new_var[valid_mask]
    filtered_adata = adata_copy[:, filtered_var.index].copy()
    filtered_adata.var = filtered_var
    
    print(f"Returning filtered AnnData with {filtered_adata.n_vars} genes")
    return filtered_adata
