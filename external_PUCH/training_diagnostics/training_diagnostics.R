#!/usr/bin/env Rscript
# Training-only diagnostics of the fixed 733-gene rank model. Uses base R and jsonlite.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: training_diagnostics.R PACKAGE_ROOT PREPARED_WORKSPACE OUTPUT_DIR")
package <- args[1]; workspace <- args[2]; output <- args[3]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
frozen <- jsonlite::read_json(file.path(package, "external_PUCH/extension_v21/results/frozen_models.json"),
                              simplifyVector = TRUE)
universe <- read.csv(file.path(package, "results/reference/universe_historical733.csv"))$gene
stopifnot(length(universe) == 733, !anyDuplicated(universe))
members <- frozen$gene_sets
ifn <- members$IFN4_rank
apm <- members$APM6_rank
hgnc <- read.delim(file.path(package, "resources/hgnc_complete_set.txt"), sep = "\t", quote = "",
                   check.names = FALSE, colClasses = "character")
hgnc <- subset(hgnc, status == "Approved" & !is.na(entrez_id) & entrez_id != "")
unique_ids <- function(x) !duplicated(x) & !duplicated(x, fromLast = TRUE)
hgnc <- hgnc[unique_ids(hgnc$symbol) & unique_ids(hgnc$entrez_id), ]
paths <- c(GSE91061 = "melanoma_revision_work/external_audit/riaz_rld.csv",
           GSE78220 = "v7_work/new_bulk/hugo2016/expression_gene_by_patient.csv",
           Liu2019 = "v7_work/extra_bulk/liu2019/expression_baseline_log2TPM1.csv.gz",
           Gide2019 = "v7_work/new_bulk/gide2019/expression_gene_by_patient.csv")
cohorts <- names(paths)
read_cohort <- function(cohort) {
    patients <- read.csv(file.path(workspace, "prepared_v7/inputs/bulk", cohort, "clinical.csv"),
                         check.names = FALSE, row.names = 1)
    patients <- patients[!is.na(patients$response), , drop = FALSE]
    expression <- as.matrix(read.csv(file.path(workspace, paths[[cohort]]), row.names = 1,
                                     check.names = FALSE))
    if (cohort == "GSE91061") rownames(expression) <- hgnc$symbol[match(rownames(expression), hgnc$entrez_id)]
    stopifnot(all(vapply(universe, function(g) sum(rownames(expression) == g, na.rm = TRUE) == 1L, logical(1))))
    expression <- expression[universe, rownames(patients), drop = FALSE]
    stopifnot(all(is.finite(expression)), all(patients$response %in% 0:1))
    ranks <- apply(expression, 2, rank, ties.method = "average")
    score <- function(genes) (colMeans(ranks[genes, , drop = FALSE]) - 1) / (length(universe) - 1) - 0.5
    data.frame(cohort = cohort, y = patients$response, IFN4 = score(ifn), APM6 = score(apm), row.names = NULL)
}
data <- do.call(rbind, lapply(cohorts, read_cohort))
sizes <- as.numeric(table(factor(data$cohort, levels = cohorts)))
responders <- as.numeric(tapply(data$y, factor(data$cohort, levels = cohorts), sum))
stopifnot(identical(sizes, c(49,26,102,72)), identical(responders, c(10,14,41,39)))
raw <- data[, c("IFN4", "APM6")]
stopifnot(max(abs(colMeans(raw) - frozen$predictor_mean)) < 1e-12,
          max(abs(vapply(raw, sd, numeric(1)) - frozen$predictor_sample_sd)) < 1e-12)
data[, c("IFN4", "APM6")] <- scale(raw, center = frozen$predictor_mean, scale = frozen$predictor_sample_sd)
auc <- function(score, y) {
    nr <- sum(y == 1); nn <- sum(y == 0)
    (sum(rank(score, ties.method = "average")[y == 1]) - nr * (nr + 1) / 2) / (nr * nn)
}
summarize <- function(d, name) {
    data.frame(cohort = name, n = nrow(d), responders = sum(d$y), auc_IFN4 = auc(d$IFN4, d$y),
               auc_APM6 = auc(d$APM6, d$y), rho_IFN4_APM6 = cor(d$IFN4, d$APM6, method = "spearman"),
               pearson_IFN4_APM6 = cor(d$IFN4, d$APM6), mean_zIFN4 = mean(d$IFN4), mean_zAPM6 = mean(d$APM6))
}
summaries <- do.call(rbind, c(list(summarize(data, "Pooled")),
                              lapply(cohorts, function(k) summarize(data[data$cohort == k, ], k))))
fit <- function(d, predictors, intercept, weighting, target_weight = nrow(d)) {
    indicator <- if (intercept == "common") matrix(1, nrow = nrow(d), dimnames = list(NULL, "common")) else {
        sapply(cohorts, function(k) as.integer(d$cohort == k))
    }
    design <- cbind(indicator, as.matrix(d[, predictors, drop = FALSE]))
    penalty <- c(rep(0, ncol(indicator)), rep(1, length(predictors)))
    w <- if (weighting == "equal_patients") rep(target_weight/nrow(d), nrow(d)) else {
        target_weight / (length(unique(d$cohort)) * as.numeric(table(d$cohort)[d$cohort]))
    }
    stopifnot(abs(sum(w) - target_weight) < 1e-10)
    objective <- function(b) {
        eta <- drop(design %*% b)
        sum(w * (pmax(eta, 0) + log1p(exp(-abs(eta))) - d$y * eta)) + sum(penalty * b^2) / 2
    }
    gradient <- function(b) drop(crossprod(design, w * (plogis(design %*% b) - d$y))) + penalty * b
    hessian <- function(b) {
        p <- drop(plogis(design %*% b))
        crossprod(design, (w * p * (1-p)) * design) + diag(penalty)
    }
    sol <- nlminb(rep(0, ncol(design)), objective, gradient = gradient, hessian = hessian,
                  control = list(rel.tol = 1e-14, x.tol = 1e-14))
    grad <- max(abs(gradient(sol$par)))
    stopifnot(grad < 1e-5)
    b <- setNames(sol$par, colnames(design))
    get_coef <- function(name) if (name %in% names(b)) unname(b[name]) else NA_real_
    result <- data.frame(intercept = intercept, weighting = weighting, model = paste(predictors, collapse = "+"),
                         n = nrow(d), responders = sum(d$y), total_weight = sum(w),
                         beta_IFN4 = get_coef("IFN4"), beta_APM6 = get_coef("APM6"),
                         intercept_common = get_coef("common"))
    for (k in cohorts) result[[paste0("intercept_", k)]] <- get_coef(k)
    result$max_abs_gradient <- grad
    result
}
models <- list("IFN4", "APM6", c("IFN4", "APM6"))
results <- list()
for (i in c("common", "cohort")) for (w in c("equal_patients", "equal_cohorts")) for (m in models) {
    results[[length(results) + 1L]] <- fit(data, m, i, w)
}
fits <- do.call(rbind, results)
reference <- subset(fits, intercept == "common" & weighting == "equal_patients")
for (name in c("IFN4_rank", "IFN4_APM6_rank")) {
    spec <- frozen$models[[name]]
    row <- reference[reference$model == if (name == "IFN4_rank") "IFN4" else "IFN4+APM6", ]
    estimates <- unlist(row[1,c("intercept_common", "beta_IFN4", if (name != "IFN4_rank") "beta_APM6")])
    stopifnot(max(abs(estimates - c(spec$intercept, spec$slopes))) < 1e-6)
}
omissions <- do.call(rbind, lapply(cohorts, function(k) {
    d <- data[data$cohort != k, ]
    cbind(omitted = k, fit(d, c("IFN4","APM6"), "common", "equal_patients", nrow(data)))
}))
write.csv(summaries, file.path(output, "cohort_rank_summaries.csv"), row.names = FALSE, na = "")
write.csv(fits, file.path(output, "training_model_coefficients.csv"), row.names = FALSE, na = "")
write.csv(omissions, file.path(output, "leave_one_cohort_out.csv"), row.names = FALSE, na = "")
