# Analyses and simulations

## `sim1`

When there is heterogeneity we can either do nothing (use all instruments), remove outliers, or adjust outliers for their detected alternative pathways.

We can assume knowledge of all outliers (pleiotropy known), or try to detect outliers (pleiotropy detected).

Therefore there are 5 methods:

- Raw, where the effects are estimated with IVW random effects model
- Outliers removed (pleiotropy known), where the SNPs known to have pleiotropic effects are removed
- Outliers removed (pleiotropy detected), where the RadialMR package is used to detect pleiotropy and then those are removed
- Outliers adjusted (pleiotropy known), where the SNPs known to have pleiotropic effects are used to find alternative pathways, and the effects are adjusted based on their associations with other traits
- Outliers adjusted (pleiotropy detected), where the SNPs detected (using RadialMR) to be outliers are used to find alternative pathways, and the effects are adjusted based on those.

Amongst these methods, how does statistical power, FDR and bias compare?

Vary:

- Number of SNPs that influence an alternative trait that in turn influences the outcome
- Number of SNPs that influence an alternative trait that in turn influences both the exposure and the outcome

Run 1000 simulations for each parameter.

To run the simulations (best to do it on bluecrystal3 node with 16 cores):

```
cd scripts
Rscript run_simulations1.r
```

To generate figures

```
cd scripts
Rscript analyse_simulations.r
```

