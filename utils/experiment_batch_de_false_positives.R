# utils/experiment_batch_de_false_positives.R
#
# Diagnostic experiment, NOT a test: how often does each proteomics DE mode call
# a protein significant when nothing differs between conditions?
#
#   legacy_corrected : ComBat (condition protected in its model), then limma
#                      ~ 0 + group on the corrected values
#   covariate        : limma ~ 0 + group + batch on the uncorrected values
#
# Under the null the fraction of raw p-values below alpha should be close to
# alpha. Fitting on ComBat output treats the batch variation it removed as free,
# so with conditions unevenly spread over batches its rate is expected to sit
# above alpha (Nygaard, Rodland & Hovig, Biostatistics 2016). The covariate fit
# is the reference. Results vary with the settings and the seeds, which is why
# this stays out of the test suite.
#
# Data are synthetic. Run from the repo root:
#   Rscript utils/experiment_batch_de_false_positives.R
# Needs limma and sva.

alpha     <- 0.05
n_prot    <- 2000L
n_reps    <- 20L
seeds     <- 1000L + seq_len(n_reps)
designs <- list(
    balanced   = list(B1 = c(Ctrl = 3, Trt = 3), B2 = c(Ctrl = 3, Trt = 3)),
    unbalanced = list(B1 = c(Ctrl = 5, Trt = 1), B2 = c(Ctrl = 1, Trt = 5))
)

#' Build one null data set for a group x batch layout
#'
#' @param layout Named list: one element per batch, a named count per condition.
#' @param seed Integer seed for the simulation.
#' @return A list with `expr` (proteins x samples) and `meta` (Group, Batch).
simulate_null <- function(layout, seed) {
    meta <- do.call(rbind, lapply(names(layout), function(b) {
        cnt <- layout[[b]]
        data.frame(Group = rep(names(cnt), times = cnt), Batch = b,
                   stringsAsFactors = FALSE)
    }))
    expr <- withr::with_seed(seed, {
        base  <- rnorm(n_prot, 20, 1)
        shift <- rnorm(n_prot, 1.5, 0.5)       # per-protein batch effect
        m <- matrix(rnorm(n_prot * nrow(meta), sd = 0.4), n_prot) + base
        m[, meta$Batch == "B2"] <- m[, meta$Batch == "B2"] + shift
        m
    })
    list(expr = expr, meta = meta)
}

#' Fraction of proteins with raw p < alpha for one fit
#'
#' @param expr Expression matrix.
#' @param design limma design with columns Ctrl and Trt.
#' @return A scalar proportion.
fp_rate <- function(expr, design) {
    cm  <- limma::makeContrasts(Trt - Ctrl, levels = design)
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), cm))
    mean(fit$p.value[, 1] < alpha)
}

results <- do.call(rbind, lapply(names(designs), function(dn) {
    do.call(rbind, lapply(seeds, function(s) {
        d <- simulate_null(designs[[dn]], s)
        g <- factor(d$meta$Group)
        b <- factor(d$meta$Batch)

        X_group <- stats::model.matrix(~ 0 + g)
        colnames(X_group) <- levels(g)
        X_batch <- cbind(X_group, batch_B2 = as.numeric(b == "B2"))

        corrected <- sva::ComBat(d$expr, batch = b, mod = stats::model.matrix(~ g),
                                 par.prior = TRUE)
        data.frame(design = dn, seed = s,
                   legacy_corrected = fp_rate(corrected, X_group),
                   covariate        = fp_rate(d$expr, X_batch))
    }))
}))

summary_tbl <- aggregate(cbind(legacy_corrected, covariate) ~ design, data = results,
                         FUN = function(x) round(mean(x), 4))
cat(sprintf("Proportion of null proteins with p < %.2f (mean over %d seeds; target %.2f):\n",
            alpha, n_reps, alpha))
print(summary_tbl, row.names = FALSE)
