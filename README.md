# CNAsearch

CNAsearch is a Python package for detecting Copy Number Alterations (CNAs) in single-cell RNA-seq data, particularly optimized for Pluripotent Stem Cells (PSCs).

- Computational Stem Cell Biology Final Project Spring 2025
- Nikhil Choudhary, Dhruv Dubey, Elizabeth Zuerbilis
- Johns Hopkins University EN.580.447 


## Installation

```bash
pip install cnasearch
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/yourteam/cnasearch.git
```

## Features

- Gene-based windowing for genomic regions
- Reference cell selection for baseline expression
- HMM-based CNA calling with probabilistic outputs
- Multi-window filtering for robust CNA detection
- Evaluation and threshold calibration tools

## Quick Start

```python
import scanpy as sc
import cnasearch as cna

# Load your data
adata = sc.read_h5ad("your_data.h5ad")

# Preprocess the data
adata = cna.preprocess(adata)

# Run the CNA detection pipeline
cna_summary = cna.run_cna_pipeline_robust(
    adata,
    window_size=1_000_000,
    min_genes_per_window=10,
    ref_method='variance',
    ref_frac=0.15,
    min_run=2,
    prob_thresh=0.5
)

# Add annotations to the AnnData object
adata = cna.add_cna_annotations_to_obs(
    adata,
    pd.DataFrame(cna_summary),
    window_size=1_000_000,
    min_genes_per_window=10
)

# Now adata.obs contains CNA annotations
```

## Advanced Usage

For more advanced usage, including threshold calibration:

```python
best_t = cna.calibrate_threshold(
    adata,
    gt_col='known_cnvs', # column with ground truth labels
    window_size=1_000_000,
    min_genes_per_window=10,
    ref_method='variance',
    ref_frac=0.15,
    min_run=2,
    overlap_fraction=0.5,
    threshs=np.linspace(0.3, 0.7, 9)
)
```

See the examples directory for more detailed use cases.

## Citation

If you use CNAsearch in your research, please cite:

```
Your Citation Here
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
