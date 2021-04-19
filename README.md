## Overview

* [code](https://github.com/xmuyulab/ImmuneCellsSig_Comment/code)
    - comment_NC_ICT.R: describing step-by-step analysis of the ImmuneCells.Sig performance on 4 melanoma ICT responses dataset. 
    - plotROC.ipynb: generating ROC figures. 
* [data](https://github.com/xmuyulab/ImmuneCellsSig_Comment/data)
    - ImmuneCell.sig signature: ImSig.rds
    - Four RNA-seq expression dataset: GSE78220_expressionMatrix.rds, BMS038.Pre.CountTable.normalized.log.rds, PRJEB23709_Pre_73samples_phenoData.rds, NatMed_103samples_GE_matrix.rds
    - Patients ICT response information: GSE78220_PhenoInfo2.rds, BMS038_phenoData.rds, PRJEB23709_Pre_73samples_phenoData.rds, NatMed_103samples_pData.rds
    
## Requirements
R version: R 3.6.3.
Dependencies: rocc (1.3), ggplot2 (3.3.2), viridis (0.5.1), Biobase (2.46.0), edgeR (3.28.1), limma (3.42.2), biomaRt (2.42.1), dplyr (1.0.0), cancerclass (1.30.0)

Python version: Python 3 (3.7.6). 
Dependencies: matplotlib (3.3.1), pandas (1.0.1), numpy (1.19.1), sklearn (0.22.1), cmocean (2.0)
