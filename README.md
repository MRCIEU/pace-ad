# PACE eczema analysis

This repo details the analysis aiming to assess the association between childhood eczema and cord blood DNA methylation. The project is part of the Pregnancy and Childhood Epigenetics (PACE) consortium. It involves meta-analysing EWAS associations across several PACE cohorts and then running supplementary analyses to explore the results. This repo is an update from Thomas Battrams's repo: https://github.com/thomasbattram/pace-ad

To run the analyses, follow the pipeline laid out in the folders seen below. READMEs within each folder should help explain something about the analyses. 

- [00_prerequisites](00_prerequisites)
- [01_aries](01_aries)
- [02_goya](02_goya)
- [03_qc](03_qc)
- [04_meta-analysis](04_meta-analysis)
	+ [04-1_meta-analysis-test](04-1_meta-analysis-test)
	+ [04-2_meta-analysis](04-2_meta-analysis-test2)
	+ [04-3_meta-analysis](04-3_meta-analysis-test3)
	+ [04-4_meta-analysis](04-4_meta-analysis-test)
- [05_power](05_power)
- [06_prs](06_prs)
- [07_candidate-genes](07_candidate-genes)
- [08_replication](08_replication)
- [09_eos-prs](09_eos-prs)

## Overview of pipeline

The main analysis in this project was the meta-analysis of eczema EWAS. Most of the data for this comes from external collaborators that sent their data in, but there are two cohorts for which we conducted the EWAS in-house - ARIES and GOYA. The EWAS in these cohorts make up the first two steps of the pipeline. The QC of EWAS data and meta-analysis setup comes next, followed by the meta-analysis itself along with some sensitivity analyses for the meta-analysis. As few results were found, power analyses were conducted to see what we could potentially detect. Then polygenic risk scores (PRS) were generated and fitted as covariates in a repeat of the EWAS in ARIES. Similarly, we adjusted for mQTLs in a repeat of the EWAS in ARIES. Finally we conducted an association analysis of AD with a polygenic risk score for eosinophil counts using a linear model in ALSPAC. The meta-analysis results around candidate genes, as identified by GWAS follow-up analyses by [Budu-Aggrey A et al.](https://doi.org/10.1038/s41467-023-41180-2), were compared to those in the rest of the genome. Then we attempted to replicate findings from a previous EWAS of eczema [Quraishi et al.](https://doi.org/10.1186/s13148-015-0108-y). 

## Summary statistics 

The summary statistics for the meta-analysis can be found here: https://doi.org/10.5281/zenodo.21459028.
