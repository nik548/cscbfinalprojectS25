# CNAsearch

CNAsearch is a Python package for detecting Copy Number Alterations (CNAs) in single-cell RNA-seq data, particularly optimized for Pluripotent Stem Cells (PSCs).

- Computational Stem Cell Biology Final Project Spring 2025
- Nikhil Choudhary, Dhruv Dubey, Elizabeth Zuerbilis
- Johns Hopkins University EN.580.447


## Installation

```bash
!pip install git+https://github.com/nik548/cscbfinalprojectS25.git
import cnasearch as cna
```

## Features

- Genomic windows constructed across chromosomes
- Reference cell selection for baseline expression
- HMM-based CNA calling with probabilistic outputs
- Multi-window filtering for robust CNA detection
- Evaluation and threshold calibration tools

## Quick Start

### Setup
```python
!pip install scanpy python-igraph leidenalg scipy umap-learn anndata hmmlearn

from google.colab import drive
drive.mount('/content/drive')
%cd /content/drive/My Drive/CSCB_Final/
%ls

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from sklearn.metrics import f1_score
import seaborn as sns
from typing import List, Dict, Optional, Tuple
from scipy.stats import zscore
from scipy.sparse import csr_matrix, issparse
import re
from sklearn.cluster import KMeans
from hmmlearn import hmm

!pip install git+https://github.com/nik548/cscbfinalprojectS25.git
import cnasearch as cna
```
### Apply Pipeline to Detect CNAs
```python
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

# Now adata.obs contains CNA annotations ("cna_regions")
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

If you use CNAsearch in your research, please cite us.

## License

This project is constructed by undergraduate students at Johns Hopkins University.
