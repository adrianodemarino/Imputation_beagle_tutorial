# Imputation beagle tutorial
Throughout the protocol we assume Bash shell.

This is a tutorial forked from: https://www.protocols.io/run/genotype-imputation-workflow-v3-0-xbgfijw

Split tutorial step by step.

1. Installation anaconda
2. Installation software
3. 

## Create new env that we call imputation

`conda create -n imputation python=3.6 anaconda`

Activate your new env:

`conda activate imputation`

Installation of the required pack/software:
```
conda install -c bioconda eagle
conda install -c bioconda beagle
conda install -c r r-base
conda install -c bioconda bcftools
conda install -c conda-forge r-data.table
conda install -c conda-forge r-sm
```

