# genecnv

Adaptive gene‑centric CNV calling pipeline for single‑cell data.

- Computational Stem Cell Biology Final Project Spring 2025
- Nikhil Choudhary, Dhruv Dubey, Elizabeth Zuerbilis
- Johns Hopkins University EN.580.447
- 
## Overview

**genecnv** provides a complete pipeline to detect copy number variations (CNVs) at the gene level in single-cell RNA‑seq datasets. It implements a novel _adaptive gene‑centric binning_ strategy that partitions the genome into bins of approximately equal gene count, preserving local genomic context while enabling robust statistical smoothing and HMM segmentation.

### Key Features

- **Adaptive gene‑centric binning**: dynamically groups genes into sequential bins ordered by chromosome and position, ensuring balanced signal aggregation and preserving genomic structure.  
- **Expression binning**: computes per-cell average expression across adaptive bins.  
- **Log‑decay smoothing**: weights neighboring bins by genomic distance using a logarithmic decay kernel.  
- **Reference centering**: centers smoothed signals against a user‑defined reference cell population or fraction.  
- **HMM segmentation & run‑length filtering**: calls gains/losses with a 3‑state Gaussian HMM and removes isolated calls below a minimum run length.  
- **Gene annotation**: integrates with MyGene.info to annotate gene coordinates when not provided.  
- **End‑to‑end pipeline**: single function to execute preprocessing through CNV call annotation.  

## Installation

This package is not on PyPI. Install directly from GitHub:

```bash
pip install git+https://github.com/nik548/cscbfinalprojectS25.git

Then
```python
import genecnv as genecnv

## Quick Start

### Setup
```python
!pip install scanpy python-igraph leidenalg scipy umap-learn anndata hmmlearn mygene

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
import mygene

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
