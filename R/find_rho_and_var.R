#' Compute pairwise correlations with optional robust / rank methods
#' @noRd
compute_cor_matrix <- function(x, method = c("spearman", "pearson", "bicor"),
                               nThreads = 1L) {
    method <- match.arg(method)
    if (method == "spearman") {
        ranks <- matrixStats::colRanks(x, ties.method = "average",
                                       preserveShape = TRUE)
        return(WGCNA::cor(ranks, nThreads = nThreads))
    }
    if (method == "pearson") {
        return(WGCNA::cor(x, nThreads = nThreads))
    }
    # bicor: biweight midcorrelation (robust)
    WGCNA::bicor(x, nThreads = nThreads, quick = 0, pearsonFallback = "none")
}

log_progress <- function(msg) {
    message(glue::glue("{date()} => {msg}"))
}

#' @title Run bootstrapping of correlations within a dataset
#' @description This function provides the more low-level functionality
#' of bootstrapping correlations of the columns within a dataset.
#' Only use this function if you want
#' a low-level interface, else \code{\link{run_csd}}
#' provides a more streamlined approach if you want to do a CSD analysis.
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom glue glue
#' @importFrom WGCNA cor bicor
#' @importFrom matrixStats colRanks
#' @param x Numeric matrix, the gene expression matrix to analyse.
#' Genes are in columns, samples are in rows.
#' @param method Correlation method: \code{"spearman"} (default),
#' \code{"pearson"}, or \code{"bicor"} (biweight midcorrelation).
#' @param shrink Logical or numeric. If \code{TRUE}, apply Schäfer–Strimmer
#' shrinkage of the mean correlation matrix toward the identity using the
#' bootstrap variances. If a number in \eqn{[0,1]}, that intensity is used.
#' If \code{FALSE} (default), no shrinkage is applied.
#' @param covariates Optional sample-level covariates to residualize out
#' before correlation (see \code{\link{residualize_expression}}).
#' @inheritParams run_csd
#' @return A list with fields \describe{
#' \item{rho}{Numeric matrix containing the bootstrapped
#'  mean of the correlation between each column}
#' \item{var}{Numeric matrix containing the bootstrapped
#'  variance of the correlation between each column}
#' \item{shrink_intensity}{The shrinkage intensity applied (0 if none)}
#' }
#' @examples
#' data("normal_expression")
#' cor_res <- run_cor_bootstrap(
#'     x = normal_expression,
#'     n_it = 100, nThreads = 2L
#' )
#' @export
run_cor_bootstrap <- function(x, n_it = 20L, nThreads = 1L,
                              verbose = TRUE, iterations_gap = 1L,
                              method = c("spearman", "pearson", "bicor"),
                              shrink = FALSE,
                              covariates = NULL) {
    method <- match.arg(method)
    blas_set_num_threads(nThreads)
    if (!is.null(covariates)) {
        if (verbose) {
            log_progress("Residualizing expression against covariates...")
        }
        x <- residualize_expression(x, covariates)
    }
    n_res <- ncol(x)
    gene_names <- colnames(x)
    rho_matrix <- array(0, c(n_res, n_res),
        dimnames = list(gene_names, gene_names)
    )
    var_matrix <- array(0, c(n_res, n_res),
        dimnames = list(gene_names, gene_names)
    )
    for (i in seq_len(n_it)) {
        if (verbose && i %% iterations_gap == 0) {
            log_progress(glue("Running bootstrap iteration {i} of {n_it}..."))
        }
        bootstrap_ind <- sample(seq_len(nrow(x)), replace = TRUE)
        cor_matrix <- compute_cor_matrix(x[bootstrap_ind, , drop = FALSE],
                                         method = method, nThreads = nThreads)
        welford_update(
            mu = rho_matrix, var = var_matrix,
            cor_matrix = cor_matrix,
            iteration = i, n_threads = nThreads
        )
    }
    var_matrix <- var_matrix / (n_it - 1)
    shrink_intensity <- 0
    if (!identical(shrink, FALSE)) {
        intensity <- if (isTRUE(shrink)) NULL else as.numeric(shrink)
        shrunk <- shrink_correlation(rho_matrix, var_rho = var_matrix,
                                     intensity = intensity)
        rho_matrix <- shrunk$rho
        shrink_intensity <- shrunk$intensity
        if (verbose) {
            log_progress(glue(
                "Applied correlation shrinkage with intensity {round(shrink_intensity, 4)}"
            ))
        }
    }
    list(rho = rho_matrix, var = var_matrix,
         shrink_intensity = shrink_intensity, method = method)
}

validate_csd_input <- function(x_1, x_2) {
    if (anyNA(x_1)) {
        stop("Argument x_1 has missing values")
    }
    if (anyNA(x_2)) {
        stop("Argument x_2 has missing values")
    }
    if (is.null(colnames(x_1)) || is.null(colnames(x_2))) {
        stop("The input matrices must be labelled with gene names")
    }
    genes_to_analyze <- intersect(colnames(x_1), colnames(x_2))
    shared_x_1 <- intersect(colnames(x_1), genes_to_analyze)
    if (length(shared_x_1) != ncol(x_1)) {
        warning(glue("{ncol(x_1)-length(shared_x_1)} genes
                    for the first condition were not found in the second"))
    }
    shared_x_2 <- intersect(colnames(x_2), genes_to_analyze)
    if (length(shared_x_2) != ncol(x_2)) {
        warning(glue("{ncol(x_2)-length(shared_x_2)} genes
                    for the second condition were not found in the first"))
    }
    list(x_1 = x_1[, genes_to_analyze, drop = FALSE],
         x_2 = x_2[, genes_to_analyze, drop = FALSE],
         genes_to_analyze = genes_to_analyze)
}

#' @title Run CSD analysis
#' @description This function implements the CSD algorithm based on the
#' one presented by Voigt et al. 2017, extended with optional robust
#' correlation measures, confounder residualization, and correlation
#' shrinkage. All pairs of genes are
#' first compared within each condition by the chosen correlation
#' and the correlation and its variance are
#' estimated by bootstrapping.
#' Finally, the results for the two conditions are
#' compared and C-, S- and D-values are computed
#' and returned.
#' @param x_1 Numeric matrix, the gene expression matrix for
#' the first condition.
#' Genes are in columns, samples are in rows.
#' The columns must be named with the name of the genes.
#' Missing values are not allowed.
#' @param x_2 Numeric matrix, the gene expression matrix for
#' the second condition.
#' @param n_it Integer, number of bootstrap iterations
#' @param nThreads Integer, number of threads to use for computations
#' @param verbose Logical, should progress be printed?
#' @param iterations_gap If output is verbose - Number of iterations between
#' each status message
#' (Default=1 - Displayed only if \code{verbose=TRUE})
#' @param method Correlation method: \code{"spearman"} (default),
#' \code{"pearson"}, or \code{"bicor"}.
#' @param shrink Logical or numeric shrinkage intensity for the mean
#' correlation matrices (see \code{\link{shrink_correlation}}).
#' @param covariates_1,covariates_2 Optional covariates for residualization
#' within each condition (same number of rows as the corresponding
#' expression matrix).
#' @return A \code{data.frame} with
#' the additional class attribute \code{csd_res} with the
#'  results of the CSD analysis.
#'  This frame has a row for each pair of genes and has the
#'  following columns: \describe{
#'  \item{\code{Gene1}}{Character, the name of the first gene}
#'  \item{\code{Gene2}}{Character, the name of the second gene}
#'  \item{\code{rho1}}{Mean correlation of the two genes in the first condition}
#'  \item{\code{rho2}}{Mean correlation of the two genes in the second condition}
#'  \item{\code{var1}}{The estimated variance of \code{rho1}
#'   determined by bootstrapping}
#'  \item{\code{var2}}{The estimated variance of \code{rho2}
#'   determined by bootstrapping}
#'  \item{\code{cVal}}{Numeric, the conserved score.
#'   A high value indicates that the co-expression of the two genes
#'   have the same sign in both conditions}
#'   \item{\code{sVal}}{Numeric, the specific score.
#'   A high value indicates that the co-expression of the two genes
#'   have a high degree of co-expression in
#'   one condition, but not the other.}
#'   \item{\code{dVal}}{Numeric, the differentiated score.
#'   A high value indicates that the co-expression of the two genes have
#'   a high degree of co-expression in
#'   both condition, but the sign of co-expression is different.}
#' }
#' @details The gene names in \code{x_1} and \code{x_2}
#' do not need to be in the same order,
#' but must be in the same namespace.
#' Only genes present in both datasets will be considered for the analysis.
#' Number of samples may differ between conditions.
#' The parallelism gained by `nThreads` applies to the computations
#' within a single iteration. The iterations are run is serial in order
#' to reduce the memory footprint.
#' @references Voigt A, Nowick K and Almaas E
#' 'A composite network of conserved and tissue specific
#' gene interactions reveals possible genetic interactions in glioma'
#' In: \emph{PLOS Computational Biology} 13(9): e1005739.
#' (doi: \url{https://doi.org/10.1371/journal.pcbi.1005739})
#' @examples
#' data("sick_expression")
#' data("normal_expression")
#' cor_res <- run_csd(
#'     x_1 = sick_expression, x_2 = normal_expression,
#'     n_it = 100, nThreads = 2L
#' )
#' c_max <- max(cor_res$cVal)
#' @export
run_csd <- function(x_1, x_2, n_it = 20L, nThreads = 1L,
                    verbose = TRUE, iterations_gap = 1L,
                    method = c("spearman", "pearson", "bicor"),
                    shrink = FALSE,
                    covariates_1 = NULL,
                    covariates_2 = NULL) {
    method <- match.arg(method)
    validated_matrices <- validate_csd_input(x_1, x_2)
    x_1 <- validated_matrices$x_1
    x_2 <- validated_matrices$x_2
    genes_to_analyze <- validated_matrices$genes_to_analyze
    if (verbose) {
        log_progress(glue("Running CSD with
\t                  {nrow(x_1)} samples from condition 1
\t                  {nrow(x_2)} samples from condition 2
\t                  Number of genes: {length(genes_to_analyze)}
\t                  Number of bootstrap iterations: {n_it}
\t                  Correlation method: {method}
\t                  Number of threads: {nThreads}"))
    }
    if (verbose) {
        log_progress("Running correlation bootstrapping on first condition...")
    }
    res_1 <- run_cor_bootstrap(x_1,
        n_it = n_it,
        nThreads = nThreads, verbose = verbose,
        iterations_gap = iterations_gap,
        method = method, shrink = shrink,
        covariates = covariates_1
    )
    if (verbose) {
        log_progress("Running correlation bootstrapping on second condition...")
    }
    res_2 <- run_cor_bootstrap(x_2,
        n_it = n_it,
        nThreads = nThreads, verbose = verbose,
        iterations_gap = iterations_gap,
        method = method, shrink = shrink,
        covariates = covariates_2
    )
    if (verbose) {
        log_progress("Summarizing results")
    }
    csd_df <- summarizeResults(
        res_1 = res_1, res_2 = res_2,
        n_threads = nThreads
    )
    # Approximate two-sided normal p-values treating C/S/D as |Z|
    csd_df$cP <- 2 * stats::pnorm(-csd_df$cVal)
    csd_df$sP <- 2 * stats::pnorm(-csd_df$sVal)
    csd_df$dP <- 2 * stats::pnorm(-csd_df$dVal)
    csd_df$cFDR <- stats::p.adjust(csd_df$cP, method = "BH")
    csd_df$sFDR <- stats::p.adjust(csd_df$sP, method = "BH")
    csd_df$dFDR <- stats::p.adjust(csd_df$dP, method = "BH")
    attr(csd_df, "method") <- method
    attr(csd_df, "shrink_intensity") <- c(
        condition_1 = res_1$shrink_intensity,
        condition_2 = res_2$shrink_intensity
    )
    class(csd_df) <- c("csd_res", class(csd_df))
    if (verbose) {
        log_progress("Returning from CSD procedure...")
    }
    csd_df
}
