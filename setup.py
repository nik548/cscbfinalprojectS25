from setuptools import setup, find_packages

setup(
    name="cnasearch",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scanpy",
        "anndata",
        "scipy",
        "scikit-learn",
        "hmmlearn",
        "matplotlib",
        "seaborn",
        "python-igraph",
        "leidenalg",
        "umap-learn",
        "mygene",  # For gene annotation
    ],
    author="Your Team Name",
    author_email="your.email@example.com",
    description="CNAsearch: A method for detecting Copy Number Alterations in single-cell data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourteam/cnasearch",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
)
