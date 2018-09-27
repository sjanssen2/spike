# Spike
Spike is one of [SNuPy's](https://snupy-aqua.bio.inf.h-brs.de/) brothers.

NGS pipeline for SNV detection in tumor context

## Install
  1. make sure you have installed Miniconda: https://conda.io/miniconda.html
  2. clone this repo into suitable location.
  3. install snakemake `conda install -c bioconda -c conda-forge snakemake`
    - if within UKD network, you might have to force conda to **not** use certificates `conda config --set ssl_verify false `
  4. install seaborn: `conda install seaborn`
