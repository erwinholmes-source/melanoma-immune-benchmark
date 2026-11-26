#!/usr/bin/env Rscript
# Independent PUCH rank-feature reconstruction and U statistic for the post hoc contrast.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: verify_external_rank.R PACKAGE_ROOT PUCH_INPUT_DIR RESULTS_DIR")
package <- args[1]; inputs <- args[2]; results <- args[3]
frozen <- jsonlite::read_json(file.path(package, "external_PUCH/extension_v21/results/frozen_models.json"),
                             simplifyVector = TRUE)
universe <- read.csv(file.path(package, "results/reference/universe_historical733.csv"))$gene
expr <- as.matrix(read.csv(file.path(inputs, "mel_puch_exp_data.csv"), row.names = 1, check.names = FALSE))
stopifnot(length(universe) == 733, identical(dim(expr), c(20293L,55L)), !anyDuplicated(rownames(expr)))
ranks <- apply(expr[universe, ,drop = FALSE], 2, rank, ties.method = "average")
make_score <- function(genes) (colMeans(ranks[genes, , drop = FALSE]) - 1) / (length(universe) - 1) - 0.5
x <- cbind(make_score(frozen$gene_sets$IFN4_rank), make_score(frozen$gene_sets$APM6_rank))
x <- scale(x, center = frozen$predictor_mean, scale = frozen$predictor_sample_sd)
fits <- read.csv(file.path(results, "training_model_coefficients.csv"))
fits <- subset(fits, intercept == "cohort" & weighting == "equal_patients")
single <- fits[fits$model == "IFN4",]
joint <- fits[fits$model == "IFN4+APM6",]
stopifnot(nrow(single) == 1L, nrow(joint) == 1L)
single_rank <- rank(x[,1] * single$beta_IFN4, ties.method = "average")
joint_rank <- rank(drop(x %*% c(joint$beta_IFN4, joint$beta_APM6)), ties.method = "average")
clin <- read.csv(file.path(inputs, "mel_puch_cli_data.csv"), row.names = 1, check.names = FALSE)
y <- clin[colnames(expr), "response"] == 1L
stopifnot(sum(y) == 14L, sum(!y) == 41L)
auc <- function(score_ranks) (sum(score_ranks[y]) - sum(y)*(sum(y)+1)/2) / (sum(y)*sum(!y))
est <- c(auc(single_rank), auc(joint_rank), auc(joint_rank)-auc(single_rank))
reported <- read.csv(file.path(results, "external_cohort_intercept_rank_auc.csv"))
stopifnot(identical(reported$model_or_contrast, c("IFN4","IFN4+APM6","joint_minus_IFN4")))
err <- max(abs(est - reported$auc))
stopifnot(err < 1e-12)
report <- list(status = "passed", n = length(y), responders = sum(y), max_abs_point_error = err,
               verified = "Both AUCs and their paired point difference from external features rebuilt in R.",
               R_version = R.version.string)
jsonlite::write_json(report, file.path(results, "external_rank_R_verification.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = 15)
