# FAHMR: Food Allergy & Health Mendelian Randomization Workflow
FAHMR is an integrated analytical pipeline that combines the High-Definition Likelihood (HDL) framework with two-sample Mendelian randomization (TSMR). Leveraging large-scale GWAS summary statistics, this pipeline was developed to investigate the genetic correlations between polygenically predicted food allergy, especially shrimp allergy and a broad spectrum of common human diseases in the East Asian populations. Moreover, FAHMR enables the assessment of the genetically predicted causal effects of food allergies on multiple disease outcomes.
## Overview
### GWAS Summary Statistics Source
GWAS summary statistics for self-reported shrimp allergy were obtained from a previously published study in 11,011 Japanese women (Khor et al. 2018. Sci Rep. doi:10.1038/s41598-017-18241-w). GWAS summary statistics for generalized food allergy (FA) and human common diseases were obtained from NBDC Human Database (Dataset ID: hum0197.v3.gwas.v1).
### GWAS Summary Statistics Preprocessing
The input GWAS summary statistics used in this pipeline must be preprocessed as follows:\
(1) Excluding SNPs located on sex chromosomes.\
(2) Removing duplicated SNPs.\
(3) Standardizing column names to match the required format: chr, pos, rsid, A1 (effect allele), A2 (reference allele), beta, se, p, eaf, N.\
An example of properly formatted GWAS summary statistics is provided in the repository for reference.\

**FAHMR comprises two main phases:**
### 1.	Genetic Correlation Evaluation by HDL Program
We integrated the collective effects of over 1.2 million representative human genome reference SNPs and employed the High-Definition Likelihood (HDL) program to assess molecular genetic correlations between FA and multiple organ system diseases.
### 2.	Causal Inference by TSMR
The primary MR analyses were conducted in two sequential steps: \
(1) First, we treated shrimp allergy, specific to the East Asian population, as the exposure and various diseases across human organ systems as the outcomes. Subsequently, the univariable MR analyses were performed with instrument variables (IVs) selection criteria: r<sup>2</sup> < 0.1, distance = 1000 kb, _P_ < 5×10<sup>-6</sup>. These analyses aimed to assess whether genetic liability to SA would affect the risk of diseases in human organ systems. \
(2) In the second step, generalized FA, specific to the Asian population, was considered as the exposure while human organ systems diseases as the outcomes to evaluate whether generalized FA would genetically affect the incidence of diseases across human organ systems. Similarly, we performed univariable MR analyses with the same thresholds mentioned in the first step. \
FAHMR incorporates rigorous MR method selection and comprehensive sensitivity analyses (e.g., heterogeneity tests, pleiotropy assessments). Full details of the MR framework and model justification can be found in the accompanying publication.
