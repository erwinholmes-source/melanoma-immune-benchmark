# Recompute selected random-list AUCs directly from baseline matrices, using base R.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript verify_draws.R prepared_v7 results_dir")
prepared <- args[1]
out <- args[2]
checks <- read.csv(file.path(out, "spot_check_sets.csv"), check.names = FALSE)
stopifnot(nrow(checks) == 40L)
cohorts <- c("GSE91061", "GSE78220", "Liu2019", "Gide2019")
replayed <- lapply(cohorts, function(cohort) {
  clinical <- read.csv(file.path(prepared, "inputs", "bulk", cohort, "clinical.csv"),
                       row.names = 1, check.names = FALSE)
  expression <- read.csv(file.path(prepared, "inputs", "rank", cohort,
                                   "source_symbol_expression.csv.gz"),
                         row.names = 1, check.names = FALSE)
  gene_lists <- strsplit(checks$genes, ";", fixed = TRUE)
  selected <- unique(unlist(gene_lists))
  expr <- as.matrix(expression[selected, rownames(clinical), drop = FALSE])
  stopifnot(all(is.finite(expr)))
  normalized <- sweep(sweep(expr, 1, rowMeans(expr), "-"), 1, apply(expr, 1, sd), "/")
  response <- clinical$response
  keep <- !is.na(response)
  response <- response[keep]
  aucs <- vapply(seq_len(nrow(checks)), function(i) {
    genes <- gene_lists[[i]]
    weight <- ifelse(genes == checks$negative_gene[i], -1, 1)
    score <- colMeans(normalized[genes, keep, drop = FALSE] * weight)
    pairs <- outer(score[response == 1], score[response == 0], "-")
    mean((pairs > 0) + 0.5 * (pairs == 0))
  }, numeric(1))
  error <- max(abs(aucs - checks[[cohort]]))
  stopifnot(is.finite(error), error < 1e-12)
  data.frame(cohort = cohort, n_draws_checked = length(aucs), max_abs_difference = error)
})
report <- do.call(rbind, replayed)
write.csv(report, file.path(out, "independent_R_checks.csv"), row.names = FALSE)
print(report)
