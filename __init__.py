"""
CNAsearch: A method for detecting Copy Number Alterations in single-cell data
"""

__version__ = "0.1.0"

# Import main functions for package-level access
from .preprocessing import preprocess
from .genomics import annotate_genes_mygene
from .windows import (
    window_genes_by_distance,
    compute_window_expression,
    select_reference_cells
)
from .hmm import (
    smooth_calls,
    call_cnas_hmm_viterbi_by_chrom,
    fit_and_get_posteriors
)
from .pipeline import (
    run_cna_pipeline_robust,
    add_cna_annotations_to_obs
)
from .evaluation import (
    calibrate_threshold,
    compare_cell_level_cnas
)
