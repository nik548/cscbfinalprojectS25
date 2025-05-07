"""
Genomic annotation utilities for genecnv.
"""

import anndata as ad
import pandas as pd
from typing import List


def annotate_genes_mygene(
    adata: ad.AnnData,
    gene_col: str = None,
    scopes: list = ['symbol', 'ensembl.gene', 'entrezgene'],
    drop_missing: bool = True
) -> ad.AnnData:
    """
    Annotate genes with genomic coordinates using MyGene.info API.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with genes to annotate
    gene_col : str, default=None
        Column in adata.var to use for gene identification. If None, uses the index.
    scopes : list, default=['symbol', 'ensembl.gene', 'entrezgene']
        Fields to search for gene IDs (see MyGene.info documentation)
    drop_missing : bool, default=True
        Whether to drop genes that lack coordinate information
    
    Returns
    -------
    AnnData
        AnnData with chromosome, start, end annotations in var
    """
    try:
        import mygene
    except ImportError:
        raise ImportError("Please install the mygene package: pip install mygene")
    
    mg = mygene.MyGeneInfo()
    
    # Determine which identifiers to use for querying
    if gene_col is None:
        # Use index as gene identifiers
        var_copy = adata.var.reset_index()
        gene_identifiers = var_copy.iloc[:, 0].tolist()
        identifier_col = var_copy.columns[0]
    else:
        # Use specified column
        if gene_col not in adata.var.columns:
            raise ValueError(f"Column '{gene_col}' not found in adata.var")
        var_copy = adata.var.reset_index()
        gene_identifiers = adata.var[gene_col].tolist()
        identifier_col = gene_col
    
    # Clean gene identifiers (remove suffixes)
    clean_identifiers = []
    original_to_clean = {}
    
    for gene_id in gene_identifiers:
        if gene_id is None or pd.isna(gene_id):
            clean_identifiers.append("")
            continue
            
        # Convert to string
        gene_id = str(gene_id)
        
        # Remove suffixes like -1, -2, -3
        import re
        clean_id = re.sub(r'-\d+$', '', gene_id)
        
        # Store mapping from original to clean
        original_to_clean[gene_id] = clean_id
        clean_identifiers.append(clean_id)
    
    # Create a working copy of var
    var = adata.var.copy()

    # --- Patch: ensure columns are object dtype to allow new categories ---
    for col in ['chromosome','start','end']:
        if col in var.columns and pd.api.types.is_categorical_dtype(var[col]):
            var[col] = var[col].astype(object)
    # End patch
    
    # Initialize coordinate columns if they don't exist
    for col in ['chromosome', 'start', 'end']:
        if col not in var.columns:
            var[col] = None

    # Find genes that need coordinates
    missing = var[['chromosome', 'start', 'end']].isnull().any(axis=1)
    genes_to_query = [clean_identifiers[i] for i, is_missing in enumerate(missing) if is_missing]
    
    # Remove empty strings or None values
    genes_to_query = [g for g in genes_to_query if g and not pd.isna(g)]
    
    if genes_to_query:
        print(f"Querying MyGene.info for {len(genes_to_query)} genes...")
        
        # Batch query to MyGene
        try:
            results = mg.querymany(
                genes_to_query,
                scopes=scopes,
                fields=['chromosome', 'genomic_pos', 'genomic_pos_hg19'],
                returnall=True
            )['out']
        except Exception as e:
            print(f"Error querying MyGene: {str(e)}")
            results = []
        
        # Process results and update var
        gene_to_coords = {}
        
        for result in results:
            query = result.get('query', '')
            if 'notfound' in result and result['notfound']:
                continue
                
            # Try both genomic_pos (GRCh38) and genomic_pos_hg19
            gp = result.get('genomic_pos')
            if gp is None:
                gp = result.get('genomic_pos_hg19')
                
            # Skip if no genomic position found
            if gp is None:
                continue
                
            # Handle different structures of genomic_pos
            if isinstance(gp, list):
                # Take the first position if it's a list
                gp = gp[0] if gp else None
                
            # Skip if genomic position is not a dictionary
            if not isinstance(gp, dict):
                continue
                
            # Extract coordinates
            chrom = result.get('chromosome') or gp.get('chr')
            start = gp.get('start')
            end = gp.get('end')
            
            # Skip if any coordinate is missing
            if None in (chrom, start, end):
                continue
                
            # Store coordinates for this gene
            gene_to_coords[query] = (chrom, start, end)
        
        # Update var with coordinates
        for i, row in var.reset_index().iterrows():
            gene_id = row[identifier_col]
            clean_id = original_to_clean.get(gene_id, gene_id)
            
            if clean_id in gene_to_coords:
                chrom, start, end = gene_to_coords[clean_id]
                var.at[gene_id, 'chromosome'] = chrom
                var.at[gene_id, 'start'] = start
                var.at[gene_id, 'end'] = end
        
        still_missing = var[['chromosome', 'start', 'end']].isnull().any(axis=1).sum()
        print(f"After query, {still_missing} genes remain un-annotated.")
    else:
        print("All genes already have genomic coordinates - skipping MyGene query.")
    
    # Filter out genes without coordinates if requested
    if drop_missing:
        var_filtered = var.dropna(subset=['chromosome', 'start', 'end'])
        print(f"Dropped {var.shape[0] - var_filtered.shape[0]} genes without coordinates.")
        
        # Create new AnnData with filtered genes
        adata_filtered = adata[:, var_filtered.index].copy()
        adata_filtered.var = var_filtered
        return adata_filtered
    else:
        # Update the original AnnData object
        adata.var = var
        return adata
