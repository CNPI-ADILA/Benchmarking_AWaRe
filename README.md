# Benchmarking AWaRe
Public repository for Benchmarking AWaRe antibiotic use code and results 

This is the repository for code and extractable results files that accompanies the manuscript: 
Cook, A. et al., Benchmarking AWaRe: estimating optimal levels of AWaRe antibiotic use in 186 countries, territories and areas based on clinical infection and resistance burden. 2026 (preprint) doi: https://doi.org/10.64898/2026.01.26.26344900

All code and datasets have been written and prepared by Aislinn Cook. For any questions please contact aicook@citystgeorges.ac.uk


# Overview of Files
The code folder includes the files used to derive antibiotic use estimates for four scenarios
1. Primary analysis
2. `High Watch` scenario
3. `Unadjusted case counts & excl. TB` scenario
4. `Alternate benchmark CTA` scenario
5. Functions code for analysis

The data folder includes the underlying datafiles that need to be read in for deriving the estimates. This includes files of the following: 
1. cluster assignments of CTAs
2. covariates for extracting total population, income and flags
3. Imputed IQVIA MIDAS dataset of AWaRe totals
4. A list of analysis countries

The results folder includes the extractable data file for each scenario as an xlsx file. 
This file is a long-format dataset with a column for AWaRe value (Total DID, Access DID, Watch DID, Reserve DID) and columns for each scenario for the median, lower and upper quantiles (95%CI) of the draws. 

