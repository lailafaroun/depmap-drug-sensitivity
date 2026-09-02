# DepMap Drug Sensitivity Prediction

## Overview
This project explores whether gene expression and mutation data from cancer cell lines can predict sensitivity to a specific drug, using public data from the [DepMap Portal](https://depmap.org/portal/download).

## Goal
Build a model that classifies cell lines as **sensitive** or **resistant** to a chosen drug, based on their molecular profile (gene expression + mutation status).

## Data Sources
All data from DepMap Portal:
- **CCLE_expression** - gene expression (TPM) per cell line
- **CCLE_mutations** - mutation status per cell line
- **Drug sensitivity** - PRISM Repurposing or GDSC (IC50 values)
- **Cell line metadata** - cancer type, tissue of origin

## Project Status
In progress - environment setup, data collection, and research

## Environment Setup
```bash
conda create -n port1 python=3.11 pandas numpy scikit-learn matplotlib seaborn scipy jupyter -y
conda activate port1
```

## Repository Structure
(to be documented as the project develops)

## Log
- **[1/September/2026]** - Project initialized, environment set up, GitHub repo created.


