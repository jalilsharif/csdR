# csdR 1.14.0

Statistical robustness extensions to the CSD pipeline:

* `run_csd()` / `run_cor_bootstrap()` gain `method` (`spearman`, `pearson`, `bicor`),
  optional Schäfer–Strimmer correlation `shrink`, and optional confounder
  residualization via `covariates_1` / `covariates_2` (or `covariates`).
* Approximate normal-tail BH FDR columns (`cFDR`, `sFDR`, `dFDR`) are attached to
  CSD results; `csd_permutation_fdr()` provides empirical condition-label
  permutation FDR when affordable.
* New helpers: `residualize_expression()`, `shrink_correlation()`,
  `annotate_csd()`, and `select_top_csd()` for labeling S/D gene pairs by name.

# csdR 1.5.1

Fixed segfault issue which did occur in `partial_argsort()` when the `n_elements`
argument was larger than the length of the input vector. In order to ensure equivalence with
`order(x, decreasing = TRUE)[1:n_elements]`, the additional elements, if any, are padded at the 
end of the answer as `NA` values.

# csdR 1.1.3

Added citation to the `csdR` article which is now printed.

# csdR 1.1.2

Made it explicit in the documentation that missing values are not allows and
wrote a test for this case.

# csdR 1.1.1

Made some minor modifications to the `README` and the vignette.

# csdR 0.99.6

* Bugfix: Open MP was not working correctly because of missing compiler flags. For this reason, the `Makevars` file has been created.
* Calculation of column ranks now uses `matrixStats::colRanks` instead of an apply statement with `base::rank`.


# csdR 0.99.0

* Added a `NEWS.md` file to track changes to the package. This is the first public version of the package.
