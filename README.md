## Fornax cluster

- 32 GPUs available at any one time.
- 5388 probes to be analysed
- 528509 SNPs
- 846 individuals

**Each probe takes 3 hours**


## Multiple testing corrections

1. Bonferroni (probes) + permutation (SNPs)
2. Correlation matrix (probes) (1 - mean(rsq) of off-diags)
3. 1000 probes permuted - maximum value
4. Number of probes required for 99th percentile of cumulative variance from PCA decomp


