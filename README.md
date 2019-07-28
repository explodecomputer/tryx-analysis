# Analyses and simulations

## Simulations

### Pleiotropy

1. The set of simulations is split across 694 chunks.

    ```
    sbatch sims_run.sh
    ```

2. Gather simulations
    
    ```
    sims_gather.sh
    ```

3. Collapse and summarise simulations

    ```
    Rscript sims_summary.r
    ```

This creates `results/sim_summary.rdata`.

4. Generate report

    ```
    Rscript -e "rmarkdown::render('sims_report.rmd', output_format='all')"
    ```


### Multivariable LASSO

- **Excluding redundant traits** - `Rscript scripts/mvmr_sims.r`



