# EI
Predict Causal Genes at GWAS Loci using Machine Learning algorithm

This directory stores codes to implement:
1. how to build summarized genemic features (required by EI) from genomic annotations matrix (provided by users).
3. how to get gene causality predictions using corresponding summarized genomic features (input of EI) and our EI model (previously trained on all 12 traits).


Files:

creatX.r: creat summarized gene features.
pred.r: perform prediction using EI models (previously trained on 12 traits).
