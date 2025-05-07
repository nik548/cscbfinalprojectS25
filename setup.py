from setuptools import setup, find_packages
from pathlib import Path

# pull in your README.md for PyPI long description
here = Path(__file__).parent
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="genecnv",
    version="0.1.0",
    author="Nikhil Choudhary",
    author_email="nchoudh5@jhu.edu",
    description="Adaptive gene-centric CNV pipeline for single-cell data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nik548/cscbfinalprojectS25",
    packages=find_packages(include=["genecnv", "genecnv.*"]),
    install_requires=[
        "numpy>=1.19",
        "pandas>=1.1",
        "scanpy>=1.9",
        "anndata>=0.7",
        "scipy>=1.5",
        "scikit-learn>=0.24",
        "hmmlearn>=0.2",
        "matplotlib>=3.3",
        "seaborn>=0.11",
        "python-igraph>=0.9",
        "mygene>=3.2",
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "genecnv=genecnv.pipeline:run_adaptive_cnv_pipeline",
        ],
    },
)
