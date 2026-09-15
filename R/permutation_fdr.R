#' @title Condition-label permutation FDR for CSD scores
#' @description Estimates empirical FDR for S- and D-scores (and optionally
#' C-scores) by pooling samples from both conditions, randomly re-assigning
#' condition labels while preserving group sizes, and recomputing CSD on each
#' permutation. This is computationally expensive; use a modest
#' \code{n_perm} and/or fewer bootstrap iterations for the null runs.
#'
#' For a cheap approximate FDR already attached to \code{\link{run_csd}}
#' results, see the \code{sFDR}/\code{dFDR} columns (normal-tail BH on the
#' C/S/D z-like scores). Those are anti-conservative when pairs are dependent;
#' this permutation procedure is more trustworthy when affordable.
#' @param x_1,x_2 Expression matrices (samples x genes), same conventions as
#' \code{\link{run_csd}}.
#' @param observed Optional pre-computed \code{run_csd} result on the true
#' labels. Computed if \code{NULL}.
#' @param n_perm Integer number of label permutations.
#' @param n_it Bootstrap iterations used for the observed result (if computed
#' here) and for each permutation.
#' @param n_it_perm Bootstrap iterations for permutations. Defaults to
#' \code{max(2L, n_it %/% 5)} for speed.
#' @param scores Character vector subset of \code{c("sVal","dVal","cVal")}
#' for which FDR is estimated.
#' @param ... Further arguments passed to \code{\link{run_csd}}
#' (e.g. \code{method}, \code{shrink}, \code{nThreads}).
#' @return The \code{observed} CSD \code{data.frame} with extra columns
#' \code{sPermFDR}, \code{dPermFDR}, and/or \code{cPermFDR}.
#' @examples
#' set.seed(1)
#' n <- 30; p <- 4
#' x1 <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("g", 1:p)))
#' x2 <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("g", 1:p)))
#' # few perms / bootstraps for a fast example only
#' csd_permutation_fdr(x1, x2, n_perm = 5L, n_it = 5L, nThreads = 1L,
#'                     verbose = FALSE)
#' @export
csd_permutation_fdr <- function(x_1, x_2, observed = NULL,
                                n_perm = 50L, n_it = 20L,
                                n_it_perm = NULL,
                                scores = c("sVal", "dVal"),
                                ...) {
    scores <- match.arg(scores, choices = c("sVal", "dVal", "cVal"),
                        several.ok = TRUE)
    if (is.null(n_it_perm)) {
        n_it_perm <- max(2L, as.integer(n_it %/% 5L))
    }
    validated <- validate_csd_input(x_1, x_2)
    x_1 <- validated$x_1
    x_2 <- validated$x_2
    genes <- validated$genes_to_analyze
    n1 <- nrow(x_1)
    n2 <- nrow(x_2)

    dots <- list(...)
    dots$verbose <- isTRUE(dots$verbose)

    if (is.null(observed)) {
        observed <- do.call(run_csd, c(
            list(x_1 = x_1, x_2 = x_2, n_it = n_it), dots
        ))
    }

    pooled <- rbind(x_1, x_2)
    null_max <- matrix(0, nrow = n_perm, ncol = length(scores),
                       dimnames = list(NULL, scores))

    for (b in seq_len(n_perm)) {
        idx <- sample.int(n1 + n2)
        p1 <- pooled[idx[seq_len(n1)], , drop = FALSE]
        p2 <- pooled[idx[n1 + seq_len(n2)], , drop = FALSE]
        colnames(p1) <- genes
        colnames(p2) <- genes
        null_res <- do.call(run_csd, c(
            list(x_1 = p1, x_2 = p2, n_it = n_it_perm), dots
        ))
        for (sc in scores) {
            null_max[b, sc] <- max(null_res[[sc]], na.rm = TRUE)
        }
    }

    # Empirical FDR via comparison to the permutation null maximum
    # (family-wise style) and a per-pair BH on empirical p-values.
    for (sc in scores) {
        obs <- observed[[sc]]
        # p = (1 + #perms with max_null >= obs_i) / (n_perm + 1)
        emp_p <- vapply(obs, function(v) {
            (1 + sum(null_max[, sc] >= v)) / (n_perm + 1)
        }, numeric(1))
        fdr <- stats::p.adjust(emp_p, method = "BH")
        p_col <- switch(sc, sVal = "sPermP", dVal = "dPermP", cVal = "cPermP")
        fdr_col <- switch(sc, sVal = "sPermFDR", dVal = "dPermFDR", cVal = "cPermFDR")
        observed[[p_col]] <- emp_p
        observed[[fdr_col]] <- fdr
    }
    attr(observed, "permutation") <- list(
        n_perm = n_perm, n_it_perm = n_it_perm,
        null_max = null_max, scores = scores
    )
    observed
}
