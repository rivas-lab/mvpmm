#!/bin/bash

# mvpmm.R arguments:
# 1. File containing summary statistics.
# 2. File to save output to.
# 3. Save to HDF5. If FALSE, save as RData file.
# 4. If mcmc=TRUE, then do an MCMC run. Otherwise, do optimization.
# 5. Number of MCMC iterations or number of draws for optimization.
# 6. Number of MCMC warmup iterations. Ignored for optimization.
# 7. Number of MCMC chains. Ignored for optimization.
# 8. Number of cores (should probably be the same as the number of chains). # Ignored for optimization.
# 9. mvpmm.stan file to use.

Rscript ../mvpmm.R \
	mvp_input.tsv.gz \
	mvp_output.RData \
	FALSE \
	FALSE \
	100 \
	0 \
	4 \
	4 \
	../mvpmm.stan
