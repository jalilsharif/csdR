#' @title Residualize expression against confounders
#' @description For each gene (column), fit a linear model on the supplied
#' covariates and replace the expression with the residuals. This is the
#' standard way to remove batch / purity / library-size effects before
#' co-expression analysis.
#' @param x Numeric matrix, samples in rows and genes in columns.
#' @param covariates Numeric matrix or data.frame with the same number of rows
#' as \code{x}. Factors are expanded via \code{model.matrix}.
#' An intercept is always included.
#' @return Numeric matrix of the same dimensions as \code{x}, with column names
#' preserved.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' colnames(x) <- paste0("g", 1:5)
#' batch <- rep(1:2, each = 5)
#' residualize_expression(x, covariates = data.frame(batch = factor(batch)))
#' @export
residualize_expression <- function(x, covariates) {
    if (is.null(covariates)) {
        return(x)
    }
    if (nrow(as.matrix(covariates)) != nrow(x)) {
        stop("covariates must have the same number of rows as x (samples)")
    }
    design <- stats::model.matrix(~ ., data = as.data.frame(covariates))
    # QR-based projection: residuals = y - Q Q^T y
    qr_design <- qr(design)
    resid_mat <- x - qr.fitted(qr_design, x)
    storage.mode(resid_mat) <- "double"
    colnames(resid_mat) <- colnames(x)
    rownames(resid_mat) <- rownames(x)
    resid_mat
}

#' @title Shrink a correlation matrix toward the identity
#' @description Applies Schäfer–Strimmer style shrinkage of off-diagonal
#' correlations toward zero. If \code{intensity} is \code{NULL}, the intensity
#' is estimated from the bootstrap variances when provided, otherwise a
#' conservative default based on the magnitude of off-diagonal entries is used.
#' @param rho Numeric square correlation matrix (mean correlations).
#' @param var_rho Optional numeric matrix of variances of the correlations
#' (e.g. from bootstrapping). Used to estimate the shrinkage intensity.
#' @param intensity Optional scalar in \eqn{[0, 1]}. If supplied, this intensity
#' is used directly instead of estimating it.
#' @return A list with \describe{
#' \item{rho}{Shrunk correlation matrix}
#' \item{intensity}{The shrinkage intensity that was applied}
#' }
#' @references Schäfer J. and Strimmer K. (2005). A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @examples
#' rho <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
#' shrink_correlation(rho, intensity = 0.2)
#' @export
shrink_correlation <- function(rho, var_rho = NULL, intensity = NULL) {
    if (!is.matrix(rho) || nrow(rho) != ncol(rho)) {
        stop("rho must be a square matrix")
    }
    n <- nrow(rho)
    if (n < 2) {
        return(list(rho = rho, intensity = 0))
    }
    off <- upper.tri(rho, diag = FALSE)
    if (is.null(intensity)) {
        if (!is.null(var_rho)) {
            num <- sum(var_rho[off])
            den <- sum(rho[off]^2)
            intensity <- if (den > 0) num / den else 1
        } else {
            # Fall back: shrink more when average |rho| is small / noisy
            mean_sq <- mean(rho[off]^2)
            intensity <- if (mean_sq > 0) min(1, 1 / (1 + n * mean_sq)) else 1
        }
        intensity <- max(0, min(1, intensity))
    } else {
        if (intensity < 0 || intensity > 1) {
            stop("intensity must be in [0, 1]")
        }
    }
    shrunk <- rho * (1 - intensity)
    diag(shrunk) <- 1
    list(rho = shrunk, intensity = intensity)
}
