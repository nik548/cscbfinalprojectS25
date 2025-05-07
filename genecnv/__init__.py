"""
genecnv: Adaptive gene‑centric CNV calling pipeline for single‑cell data.
"""

__version__ = "0.1.0"

# expose main pipeline functions
from .pipeline import run_adaptive_cnv_pipeline, annotate_cnv_calls

# expose utilities
from .binning import make_adaptive_bins
from .expression import bin_expression_matrix
from .smoothing import log_decay_smoothing
from .baseline import reference_center
from .segmentation import call_hmm, run_length_filter
from .preprocessing import preprocess
from .annotate import annotate_genes_mygene


# package‑level API
__all__ = [
    "run_adaptive_cnv_pipeline",
    "annotate_cnv_calls",
    "make_adaptive_bins",
    "bin_expression_matrix",
    "log_decay_smoothing",
    "reference_center",
    "call_hmm",
    "run_length_filter",
    "preprocess",
    "annotate_genes_mygene"
]
