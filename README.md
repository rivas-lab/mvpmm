# MultiVariate Polygenic Mixture Model (MVPMM)
## Christopher DeBoever and Manuel A. Rivas
This directory contains the MultiVariate Polygenic Mixture Model (MPVMM) for
estimating genetic parameters as described in 
[DeBoever et al. 2019](https://www.biorxiv.org/content/10.1101/738856v2).

MVPMM takes as input GWAS summary statistics (estimated effect sizes and
standard errors) for two phenotypes and estimates genetic parameters such as
genetic correlation, polygenicity, and scale of genetic effects for the
phenotypes. Please see [the manuscript](https://www.biorxiv.org/content/10.1101/738856v2) for more
information on MVPMM.
<img src="https://www.biorxiv.org/content/biorxiv/early/2019/08/21/738856/F4.large.jpg"  />

## Running MVPMM

MVPMM is specified in the `mvpmm.stan` file. We've included a helper script
`mvpmm.R` for running MVPMM using [rstan](https://cran.r-project.org/package=rstan).
`mvpmm.R` has the following arguments:

1. `data.fn`: File containing summary statistics.
data.fn
2. `out.fn`: File to save output to.
out.fn
3. `hdf5`: Save to HDF5. If FALSE, save as RData file.
4. `mcmc`: If `mcmc=TRUE`, then do an MCMC run. Otherwise, use Stan's
   optimization function to perform maximum a posteriori estimation.
5. `niter`: Number of MCMC iterations or number of draws for optimization.
6. `warmup`: Number of MCMC warmup iterations. Ignored for optimization.
7. `chains`: Number of MCMC chains. Ignored for optimization.
8. `cores`: Number of cores (should probably be the same as the number of
chains). Ignored for optimization.
9. `mvpmm.fn`: Stan file to use.

## Input File Format

An input file is a tab-separated file with an index column followed by these
columns:

* `ID`: Variant ID (e.g. rsID)
* `BETA`: GWAS effect size estimate
* `SE`: Standard error of GWAS effect size
* `P`: Unadjusted p-value from GWAS
* `code`: Phenotype label

The variant IDs contained in the `ID` column should not be repeated within a
phenotype. Each variant should be present for each phenotype. Take care to LD
prune variants prior to creating the input file and make sure that the effect
alleles are consistent between GWAS studies.

This is an example of a few lines from an input file:

```
              ID      BETA        SE         P code
1085   rs1001495 -0.072786  0.036405  0.045570   RA
6055   rs1001495  0.029559  0.051766  0.567996  SLE
1530  rs10025152 -0.021224  0.024549  0.387300   RA
6500  rs10025152 -0.020203  0.030855  0.512618  SLE
1385  rs10031922  0.075478  0.046814  0.106900   RA
6355  rs10031922  0.009950  0.103314  0.923273  SLE
```

There is an example input file included at `example/mvp_input.tsv.gz`. 

## Example

The file `example/mvp_input.tsv.gz` is an example input file created from
GWAS summary statistics for family history of high blood pressure and essential
hypertension in the UK Biobank. Note that the variants included in this file
have been LD-pruned.

We can run MVPMM using Stan's optimization function for 100 iterations for
these GWAS summary statistics using the following command:

	Rscript mvpmm.R example/mvp_input.tsv.gz mvp_output.RData FALSE FALSE 100 0 1 1 mvpmm.stan

We can perform an MCMC with 100 burn-in iterations run as follows:

	Rscript mvpmm.R example/mvp_input.tsv.gz ra_sle.RData FALSE TRUE 400 100 1 1 mvpmm.stan

You can change into the `example` directory and use `bash run.sh` to run an
optimization example.

