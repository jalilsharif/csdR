#' @title Annotate CSD results with edge labels
#' @description Assigns each gene pair a primary label based on which of the
#' C / S / D scores is largest, and for S-edges reports which condition has
#' the stronger absolute co-expression.
#' @param csd_df A \code{data.frame} as returned by \code{\link{run_csd}}.
#' @param condition_names Character vector of length 2 giving display names for
#' condition 1 and condition 2 (defaults to \code{c("condition_1", "condition_2")}).
#' @return The input frame with two extra columns: \describe{
#' \item{label}{Character, one of \code{"C"}, \code{"S"}, \code{"D"}}
#' \item{stronger_in}{For S edges, the condition with larger \eqn{|\\rho|};
#' \code{NA} otherwise.}
#' }
#' @examples
#' set.seed(1)
#' x1 <- matrix(rnorm(40), 10, 4, dimnames = list(NULL, paste0("g", 1:4)))
#' x2 <- matrix(rnorm(40), 10, 4, dimnames = list(NULL, paste0("g", 1:4)))
#' res <- run_csd(x1, x2, n_it = 5L, verbose = FALSE)
#' annotate_csd(res, condition_names = c("normal", "sick"))
#' @export
annotate_csd <- function(csd_df, condition_names = c("condition_1", "condition_2")) {
    if (!all(c("cVal", "sVal", "dVal", "rho1", "rho2") %in% names(csd_df))) {
        stop("csd_df must contain cVal, sVal, dVal, rho1 and rho2 columns")
    }
    if (length(condition_names) != 2) {
        stop("condition_names must have length 2")
    }
    score_mat <- cbind(C = csd_df$cVal, S = csd_df$sVal, D = csd_df$dVal)
    label <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
    stronger_in <- rep(NA_character_, nrow(csd_df))
    is_s <- label == "S"
    if (any(is_s)) {
        stronger_in[is_s] <- ifelse(
            abs(csd_df$rho1[is_s]) >= abs(csd_df$rho2[is_s]),
            condition_names[[1]],
            condition_names[[2]]
        )
    }
    csd_df$label <- label
    csd_df$stronger_in <- stronger_in
    csd_df
}

#' @title Select top-scoring C, S and D gene pairs
#' @description Convenience wrapper around \code{\link{partial_argsort}} that
#' returns the top pairs for each score with an explicit \code{label} column.
#' @param csd_df A \code{data.frame} as returned by \code{\link{run_csd}}.
#' @param n_pairs Integer, number of top pairs to keep for each of C, S and D.
#' @param condition_names Passed through for S-edge \code{stronger_in} labels.
#' @return A \code{data.frame} containing the union of top C/S/D pairs with a
#' \code{label} column indicating which list they were selected from. If a pair
#' appears in more than one top list, the label corresponding to its highest
#' score among those lists is kept.
#' @export
select_top_csd <- function(csd_df, n_pairs = 100L,
                           condition_names = c("condition_1", "condition_2")) {
    n_pairs <- as.integer(n_pairs)
    n_pairs <- min(n_pairs, nrow(csd_df))
    frames <- list(
        C = csd_df[partial_argsort(csd_df$cVal, n_pairs), , drop = FALSE],
        S = csd_df[partial_argsort(csd_df$sVal, n_pairs), , drop = FALSE],
        D = csd_df[partial_argsort(csd_df$dVal, n_pairs), , drop = FALSE]
    )
    for (lab in names(frames)) {
        frames[[lab]]$label <- lab
        frames[[lab]]$score_for_label <- switch(
            lab,
            C = frames[[lab]]$cVal,
            S = frames[[lab]]$sVal,
            D = frames[[lab]]$dVal
        )
    }
    out <- do.call(rbind, frames)
    key <- paste(out$Gene1, out$Gene2, sep = "\r")
    ord <- order(key, -out$score_for_label)
    out <- out[ord, , drop = FALSE]
    out <- out[!duplicated(paste(out$Gene1, out$Gene2, sep = "\r")), , drop = FALSE]
    out$score_for_label <- NULL

    out$stronger_in <- NA_character_
    is_s <- out$label == "S"
    if (any(is_s)) {
        out$stronger_in[is_s] <- ifelse(
            abs(out$rho1[is_s]) >= abs(out$rho2[is_s]),
            condition_names[[1]],
            condition_names[[2]]
        )
    }
    rownames(out) <- NULL
    out
}
