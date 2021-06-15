# TFM_analysis
Custom scripts used for the TFM thesis

This repository contains files with the code used to generate the figures present in the resuts section from the TFM thesis entitled **"Impact of NIPBL downregulation in the behaviour of cohesin variants STAG1 and STAG2: consequences for gene regulation and genome folding"**

The repository is organized according to the figure order of the written project. 

## Fig2.1: RNA seq
1. DEG analysis using `DESeq2`
2. Heatmap and Volcano Plots of gene expresion
3. GSEA 
4. Heatmap of Cohesin distibution along TSS of genes

## Fig 2.2: Super enhancers
1. Heatmap of Cohesin distribution at Super Enhancer regions
2. Odds Ratio SuperEnhancer/NoSuperEnhancer for the 3 DEG categories

## Fig 2.3: Cohesin Distribution genome wide
1. Heatmap of log2FoldChange siNIPBL/WT for SMC1, SA2 and SA1
2. Odds Ratio siNIPBL/WT of cohesin peaks at the 18 states model

## Fig 2.4: Magnitude of Cohesin change
1. Boxplot of log2FC at the simplified 6 states model
2. SA2 Score dynamics
3. SA1 Score dynamics

## Fig 2.5: Changes with respect to TAD borders and distance to loading sites
1. Boxplots log2FC with respect to TAD borders
2. log2FC with respect to BRD4 (NIPBL) peaks


