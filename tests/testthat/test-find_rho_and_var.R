set.seed(432)
n_it <- 10000
N_SAMPLES <- 100
gene_1 <- rnorm(N_SAMPLES, mean = 0, sd = 1)
gene_2 <- 1 / sqrt(2) * gene_1 + 1 / sqrt(2) * rnorm(N_SAMPLES, mean = 0, sd = 1)
combined_matrix <- cbind(gene_1, gene_2)
gene_names <- c("gene_1", "gene_2")
gene_names_dimnames <- list(gene_names, gene_names)
expected_res <- replicate(n_it, {
    bootstrap_sample <- sample(1:N_SAMPLES, size = N_SAMPLES, replace = TRUE)
    stats::cor(x = gene_1[bootstrap_sample], y = gene_2[bootstrap_sample], method = "spearman")
})
expected_rho <- mean(expected_res)
expected_var <- var(expected_res)
expected_result <- list(
    rho = matrix(c(1, expected_rho, expected_rho, 1), byrow = FALSE, ncol = 2, dimnames = gene_names_dimnames),
    var = matrix(c(0, expected_var, expected_var, 0), byrow = FALSE, ncol = 2, dimnames = gene_names_dimnames)
)
test_that("run_cor_bootstrap works as intended", {
    cor_bootstrap_res <- run_cor_bootstrap(combined_matrix, n_it = 1000, nThreads = 2L, verbose = FALSE)
    expect_equal(object = cor_bootstrap_res$rho, expected = expected_result$rho, tolerance = 0.05)
    expect_equal(object = cor_bootstrap_res$var, expected = expected_result$var, tolerance = 0.05)
})

gene_1_sick <- rnorm(N_SAMPLES, mean = 0, sd = 2)
gene_2_sick <- 1 / sqrt(2) * gene_1 - 1 / sqrt(2) * rnorm(N_SAMPLES, mean = 0, sd = 2)
combined_matrix_sick <- cbind(gene_1_sick, gene_2_sick)
colnames(combined_matrix_sick) <- gene_names
expected_res_sick <- replicate(n_it, {
    bootstrap_sample <- sample(1:N_SAMPLES, size = N_SAMPLES, replace = TRUE)
    stats::cor(x = gene_1_sick[bootstrap_sample], y = gene_2_sick[bootstrap_sample], method = "spearman")
})
expected_rho_sick <- mean(expected_res_sick)
expected_var_sick <- var(expected_res_sick)
expected_sd_estimate <- sqrt(expected_var + expected_var_sick)
expected_c <- abs(expected_rho + expected_rho_sick) / expected_sd_estimate
expected_s <- abs(abs(expected_rho) - abs(expected_rho_sick)) / expected_sd_estimate
expected_d <- abs(abs(expected_rho) + abs(expected_rho_sick) - abs(expected_rho + expected_rho_sick)) / expected_sd_estimate

test_that("Overall CSD algorithm works with a small dataset", {
    res <- run_csd(
        x_1 = combined_matrix, x_2 = combined_matrix_sick,
        n_it = n_it, nThreads = 2, verbose = FALSE
    )
    expect_equal(res$Gene1, "gene_1")
    expect_equal(res$Gene2, "gene_2")
    expect_equal(res$rho1, expected_rho, tolerance = 0.05)
    expect_equal(res$rho2, expected_rho_sick, tolerance = 0.05)
    expect_equal(res$cVal, expected_c, tolerance = 0.05)
    expect_equal(res$sVal, expected_s, tolerance = 0.05)
    expect_equal(res$dVal, expected_d, tolerance = 0.05)
    expect_true(all(c("cFDR", "sFDR", "dFDR") %in% names(res)))
})

test_that("Overall CSD algorithm gives a reasonable result for more realistic data", {
    data("normal_expression")
    data("sick_expression")
    csd_res <- run_csd(x_1 = normal_expression, x_2 = sick_expression, n_it = 20, nThreads = 2, verbose = FALSE)
    expect_equal(object = nrow(csd_res), as.integer(1000 * 999 / 2))
})

test_that("Missing values results in errors", {
    data("normal_expression")
    data("sick_expression")
    normal_expression[20L, 30L] <- NA
    expect_error(run_csd(x_1 = normal_expression, x_2 = sick_expression,
                         n_it = 20, nThreads = 2, verbose = FALSE), regexp =  "missing values")
})

test_that("pearson and bicor methods run", {
    set.seed(2)
    x1 <- matrix(rnorm(80), 20, 4, dimnames = list(NULL, paste0("g", 1:4)))
    x2 <- matrix(rnorm(60), 15, 4, dimnames = list(NULL, paste0("g", 1:4)))
    pear <- run_csd(x1, x2, n_it = 10L, method = "pearson", verbose = FALSE)
    bic <- run_csd(x1, x2, n_it = 10L, method = "bicor", verbose = FALSE)
    expect_equal(nrow(pear), 6L)
    expect_equal(nrow(bic), 6L)
    expect_equal(attr(pear, "method"), "pearson")
    expect_equal(attr(bic, "method"), "bicor")
})

test_that("shrinkage reduces off-diagonal magnitudes", {
    rho <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
    shrunk <- shrink_correlation(rho, intensity = 0.5)
    expect_equal(shrunk$intensity, 0.5)
    expect_equal(shrunk$rho[1, 2], 0.45)
    expect_equal(diag(shrunk$rho), c(1, 1))
})

test_that("residualize_expression removes linear covariate effects", {
    set.seed(3)
    n <- 40
    batch <- rep(c(0, 1), each = n / 2)
    signal <- rnorm(n)
    x <- cbind(g1 = signal + 5 * batch, g2 = rnorm(n))
    resid <- residualize_expression(x, covariates = data.frame(batch = batch))
    fit <- lm(resid[, "g1"] ~ batch)
    expect_equal(unname(coef(fit)[["batch"]]), 0, tolerance = 1e-8)
})

test_that("annotate_csd and select_top_csd label gene names", {
    set.seed(4)
    x1 <- matrix(rnorm(100), 20, 5, dimnames = list(NULL, paste0("gene", 1:5)))
    x2 <- matrix(rnorm(100), 20, 5, dimnames = list(NULL, paste0("gene", 1:5)))
    # Induce an S-like pair: correlated only in x1
    x1[, 1] <- x1[, 2] + rnorm(20, sd = 0.1)
    res <- run_csd(x1, x2, n_it = 30L, verbose = FALSE)
    ann <- annotate_csd(res, condition_names = c("ctrl", "case"))
    expect_true(all(ann$label %in% c("C", "S", "D")))
    expect_true(all(ann$Gene1 %in% colnames(x1)))
    top <- select_top_csd(res, n_pairs = 2L, condition_names = c("ctrl", "case"))
    expect_true(all(top$label %in% c("C", "S", "D")))
    expect_true(nrow(top) <= 6L)
})

test_that("csd_permutation_fdr adds permutation FDR columns", {
    set.seed(5)
    x1 <- matrix(rnorm(45), 15, 3, dimnames = list(NULL, paste0("g", 1:3)))
    x2 <- matrix(rnorm(30), 10, 3, dimnames = list(NULL, paste0("g", 1:3)))
    out <- csd_permutation_fdr(x1, x2, n_perm = 4L, n_it = 4L, n_it_perm = 2L,
                               verbose = FALSE, nThreads = 1L)
    expect_true(all(c("sPermFDR", "dPermFDR", "sPermP", "dPermP") %in% names(out)))
})
