Necromass Metabolite Profiling: Composition and Bacterial Utilization

This repository contains: 
metabolite abundance table recording normalized levels of 865 annotated metabolites across:
    Necromass composition profiling (60 samples)
    Necromass utilization profiling (65 samples)
Associated metadata for samples and metabolites

Analytical Workflow

1. Differential Abundance Analysis

To identify metabolites consumed by at least one bacterial species, the abundance of metabolites in each post-consumption necromass (spent culture) were compared with that in the twelve-species necromass. The criteria used to identify consumed metabolites were as follows: (a) FDR-adjusted p < 0.05 in the two-sample t-test, (b) variable importance in projection (VIP) score > 1.0 in orthogonal projections to latent structures-discriminant analysis (OPLS-DA), and (c) significantly lower abundance in post-consumption necromass than in the twelve-species necromass. These consumed metabolites were further subjected to K-W tests (FDR < 0.05) to examine species-specific necromass consumption. 


Output Columns:

Metabolite ID
Mean abundance (control group - fresh necromass)
Mean abundance (treatment group - spent necromass)
VIP score (OPLS-DA)
t-test p-value
FDR-adjusted p-value
Fold Change (treatment/control)
log10-transformed Fold Change


2. Kruskal-Wallis Tests
  A. Necromass Composition Analysis
     Identified significant abundance differences across 12 single-species necromass preparations (FDR < 0.05).

  B. Necromass Depletion Analysis
    Assessed interspecific variation in metabolite consumption patterns, restricted to metabolites passing the differential abundance filters above.

Implementation
All analyses were performed using R scripts included in this repository.
   
