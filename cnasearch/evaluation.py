"""
Evaluation and threshold calibration utilities for CNAsearch.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Dict, Union, Optional, List

from .pipeline import run_cna_pipeline_robust, add_cna_annotations_to_obs


def compare_cell_level_cnas(
    adata: sc.AnnData,
    gt_col: str,
    pred_col: str,
    overlap_fraction: float = 0.5
) -> Dict[str, Union[float, int]]:
    """
    Compare predicted CNAs against ground truth at cell level.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with ground truth and predictions
    gt_col : str
        Column in adata.obs with ground truth CNA annotations
    pred_col : str
        Column in adata.obs with predicted CNA annotations
    overlap_fraction : float, default=0.5
        Minimum required overlap fraction to consider a match
    
    Returns
    -------
    metrics : Dict[str, Union[float, int]]
        Dictionary with evaluation metrics:
        - precision: TP / (TP + FP)
        - recall: TP / (TP + FN)
        - specificity: TN / (TN + FP)
        - F1: 2 * (precision * recall) / (precision + recall)
        - TP, TN, FP, FN: counts of true positives, true negatives, etc.
    
    Notes
    -----
    Regions are compared based on genomic coordinates. A prediction matches
    ground truth if their genomic spans overlap by at least overlap_fraction.
    """
    # Helper function to parse region strings into sets of (chrom, start, end, type)
    def parse_regions(s):
        if not s or pd.isna(s) or s == "":
            return set()
        
        regions = set()
        for region in s.split(','):
            # Parse format like "chr1:1000000-2000000(gain)"
            loc, call = region.strip().split('(')
            call = call.rstrip(')')
            chrom, pos = loc.split(':')
            start, end = map(int, pos.split('-'))
            regions.add((chrom, start, end, call))
        
        return regions
    
    # Parse ground truth and predictions
    gt_regions = adata.obs[gt_col].fillna("").apply(parse_regions)
    pred_regions = adata.obs[pred_col].fillna("").apply(parse_regions)
    
    # Calculate TP, FP, TN, FN at cell level
    TP, FP, TN, FN = 0, 0, 0, 0
    
    for i, (gt_set, pred_set) in enumerate(zip(gt_regions, pred_regions)):
        has_gt = len(gt_set) > 0
        has_pred = len(pred_set) > 0
        
        if has_gt and has_pred:
            # Match regions by overlap
            matches = 0
            for pred_reg in pred_set:
                p_chrom, p_start, p_end, p_call = pred_reg
                
                for gt_reg in gt_set:
                    g_chrom, g_start, g_end, g_call = gt_reg
                    
                    # Check for matching chromosome and call type
                    if p_chrom != g_chrom or p_call != g_call:
                        continue
                    
                    # Calculate overlap
                    overlap_start = max(p_start, g_start)
                    overlap_end = min(p_end, g_end)
                    
                    if overlap_end < overlap_start:
                        continue  # No overlap
                    
                    overlap_len = overlap_end - overlap_start
                    pred_len = p_end - p_start
                    gt_len = g_end - g_start
                    
                    # Check if overlap fraction is sufficient
                    if (overlap_len / pred_len >= overlap_fraction and 
                        overlap_len / gt_len >= overlap_fraction):
                        matches += 1
                        break
            
            if matches > 0:
                TP += 1
            else:
                FP += 1
        elif has_pred:
            FP += 1
        elif has_gt:
            FN += 1
        else:
            TN += 1
    
    # Calculate metrics
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    F1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {
        'precision': precision,
        'recall': recall,
        'specificity': specificity,
        'F1': F1,
        'TP': TP,
        'TN': TN,
        'FP': FP,
        'FN': FN
    }


def calibrate_threshold(
    adata: sc.AnnData,
    gt_col: str,
    window_size: int,
    min_genes_per_window: int,
    ref_method: str,
    ref_frac: float,
    min_run: int,
    overlap_fraction: float = 0.5,
    threshs: np.ndarray = np.linspace(0.1, 0.9, 17)
) -> float:
    """
    Calibrate the probability threshold for optimal CNA calling.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with ground truth
    gt_col : str
        Column in adata.obs with ground truth CNA annotations
    window_size : int
        Size of each genomic window in base pairs
    min_genes_per_window : int
        Minimum number of genes required in a window
    ref_method : str
        Method for selecting reference cells
    ref_frac : float
        Fraction of cells to select as reference
    min_run : int
        Minimum number of consecutive windows required for a call
    overlap_fraction : float, default=0.5
        Minimum required overlap fraction to consider a match
    threshs : ndarray, default=np.linspace(0.1, 0.9, 17)
        Array of thresholds to test
    
    Returns
    -------
    best_thresh : float
        Threshold with highest Youden's J statistic (recall + specificity - 1)
    
    Notes
    -----
    This function runs the full pipeline multiple times with different thresholds
    and selects the threshold that maximizes Youden's J statistic.
    """
    best_t, best_J = None, -1.0
    
    print("Starting threshold calibration over values:", np.round(threshs, 3))
    for idx, t in enumerate(threshs, start=1):
        print(f"\n--- Testing threshold {t:.3f} ({idx}/{len(threshs)}) ---")
        
        # 1) Call pipeline
        print("Running pipeline...")
        summary = run_cna_pipeline_robust(
            adata,
            window_size=window_size,
            min_genes_per_window=min_genes_per_window,
            ref_method=ref_method,
            ref_frac=ref_frac,
            min_run=min_run,
            prob_thresh=t
        )
        print(f"Pipeline returned {len(summary)} total CNA events.")
        
        # 2) Annotate predictions onto a fresh copy of adata
        ad = adata.copy()
        pred_df = pd.DataFrame(summary)
        ad = add_cna_annotations_to_obs(
            ad,
            pred_df,
            window_size=window_size,
            min_genes_per_window=min_genes_per_window
        )
        ad.obs['predictions'] = ad.obs['cna_regions']
        num_pred_pos = (ad.obs['predictions'] != "").sum()
        print(f"Number of cells with ≥1 predicted event: {num_pred_pos}/{ad.shape[0]}")
        
        # 3) Compute metrics
        m = compare_cell_level_cnas(
            ad,
            gt_col=gt_col,
            pred_col='predictions',
            overlap_fraction=overlap_fraction
        )
        J = m['recall'] + m['specificity'] - 1
        print(f"Results: recall={m['recall']:.3f}, specificity={m['specificity']:.3f}, "
              f"TP={m['TP']}, TN={m['TN']}, FP={m['FP']}, FN={m['FN']}, J={J:.3f}")
        
        if J > best_J:
            best_J, best_t = J, t
            print(f"--> New best threshold: {best_t:.3f} with J={best_J:.3f}")
    
    print(f"\nBest threshold overall = {best_t:.3f} (Youden's J = {best_J:.3f})")
    return best_t
