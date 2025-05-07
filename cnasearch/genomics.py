"""
Genomic annotation utilities for CNAsearch.
"""

import anndata as ad
import pandas as pd
from typing import List


def annotate_genes_mygene(
    adata: ad.AnnData,
    species: str = 'human',
    scopes: List[str] = ['symbol', 'accession', 'ensembl.gene', 'ensembl.transcript']
) -> ad.AnnData:
    """
    Annotate genes with genomic coordinates using MyGene.info
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with gene identifiers as var_names
    species : str, default='human'
        Species for MyGene.info query
    scopes : List[str], default=['symbol', 'accession', 'ensembl.gene', 'ensembl.transcript']
        Identifier types to query
    
    Returns
    -------
    AnnData
        AnnData with chromosome, start, end annotations in var
    
    Notes
    -----
    Requires mygene package: pip install mygene
    """
    try:
        import mygene
    except ImportError:
        raise ImportError("Please install the mygene package: pip install mygene")
    
    mg = mygene.MyGeneInfo()
    
    # 1) Pull var → DataFrame, rename index to 'gene'
    idx_name = adata.var.index.name or 'index'
    var = adata.var.reset_index().rename(columns={idx_name: 'gene'})
    
    # 2) Mask genes missing ANY of chromosome/start/end
    missing = var[['chromosome', 'start', 'end']].isnull().any(axis=1)
    to_query = var.loc[missing, 'gene'].tolist()
    
    if to_query:
        print(f"Querying MyGene.info for {len(to_query)} genes…")
        # returnall=True gives {'out': [...]}
        res = mg.querymany(
            to_query,
            scopes=scopes,
            fields='chromosome,genomic_pos',
            species=species,
            returnall=True
        )['out']
        
        # 3) Build a lookup
        ann = {}
        for rec in res:
            gene = rec.get('query')
            chrom = rec.get('chromosome')
            gp = rec.get('genomic_pos')
            # unify gp into a dict
            if gp is None:
                gp_dict = {}
            elif isinstance(gp, list):
                gp_dict = gp[0] if gp else {}
            elif isinstance(gp, dict):
                gp_dict = gp
            else:
                gp_dict = {}
            
            start = gp_dict.get('start')
            end = gp_dict.get('end')
            # only save if both present
            if chrom is not None and start is not None and end is not None:
                ann[gene] = (chrom, start, end)
        
        # 4) Apply back into var DataFrame
        for i, row in var.loc[missing].iterrows():
            g = row['gene']
            if g in ann:
                var.at[i, 'chromosome'], var.at[i, 'start'], var.at[i, 'end'] = ann[g]
        
        still_missing = var[['chromosome', 'start', 'end']].isnull().any(axis=1).sum()
        print(f"After query, {still_missing} genes remain un-annotated (will be dropped).")
    else:
        print("All genes already have genomic coords—skipping MyGene query.")
    
    # 5) Drop any that still lack coords and re-assign to adata.var
    var = var.dropna(subset=['chromosome', 'start', 'end']).set_index('gene')
    adata.var = var
    return adata
