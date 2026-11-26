#!/usr/bin/env Rscript
# Independent raw-to-point-estimate check using R midranks, Mann–Whitney U and nlminb.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4)
package <- args[1]; workspace <- args[2]; puch <- args[3]; output <- args[4]
read_matrix <- function(path) as.matrix(read.csv(path, row.names=1, check.names=FALSE))
rank_auc <- function(score, y) {
    p <- sum(y==1); n <- sum(y==0)
    (sum(rank(score, ties.method="average")[y==1]) - p*(p+1)/2)/(p*n)
}
stand <- read.csv(file.path(output,"PUCH_panel_results.csv"), check.names=FALSE)
xp <- read_matrix(file.path(puch,"mel_puch_exp_data.csv"))
clin <- read.csv(file.path(puch,"mel_puch_cli_data.csv"),row.names=1,check.names=FALSE)
yext <- as.integer(clin[colnames(xp),"response"]==1)
stopifnot(length(yext)==55, sum(yext)==14)
make_z <- function(genes, panel) {
    x <- scale(t(xp[genes,,drop=FALSE]))
    if (panel=="GEP-18") x[,"CD276"] <- -x[,"CD276"]
    rowMeans(x)
}
ifn <- c("IFNG","IDO1","CXCL10","HLA-DRA")
apm <- c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAPBP")
base_score <- make_z(ifn,"IFN-4")
diff_auc <- c(); diff_rho <- c()
for(i in seq_len(nrow(stand))) {
    s <- make_z(strsplit(stand$genes_used[i],";",fixed=TRUE)[[1]], stand$panel[i])
    diff_auc[i] <- abs(rank_auc(s,yext) - stand$auc[i])
    diff_rho[i] <- abs(cor(s,base_score,method="spearman") - stand$rho_IFN4[i])
}
universe <- read.csv(file.path(package,"results/reference/universe_historical733.csv"))$gene
rank_features <- function(x) {
    ranking <- apply(x[universe,,drop=FALSE],2,rank,ties.method="average")
    means <- cbind(colMeans(ranking[ifn,]),colMeans(ranking[apm,]))
    (means-1)/(length(universe)-1)-0.5
}
hgnc <- read.delim(file.path(package,"resources/hgnc_complete_set.txt"),sep="\t",quote="",
                   check.names=FALSE,colClasses="character")
hgnc <- hgnc[hgnc$status=="Approved" & !is.na(hgnc$entrez_id) & hgnc$entrez_id!="",]
hgnc <- hgnc[!duplicated(hgnc$symbol) & !duplicated(hgnc$symbol,fromLast=TRUE) &
             !duplicated(hgnc$entrez_id) & !duplicated(hgnc$entrez_id,fromLast=TRUE),]
cohorts <- c("GSE91061","GSE78220","Liu2019","Gide2019")
files <- c("melanoma_revision_work/external_audit/riaz_rld.csv",
           "v7_work/new_bulk/hugo2016/expression_gene_by_patient.csv",
           "v7_work/extra_bulk/liu2019/expression_baseline_log2TPM1.csv.gz",
           "v7_work/new_bulk/gide2019/expression_gene_by_patient.csv")
xs <- list(); ys <- list()
for(i in seq_along(cohorts)) {
    clinical <- read.csv(file.path(workspace,"prepared_v7","inputs","bulk",cohorts[i],"clinical.csv"),
                         row.names=1,check.names=FALSE)
    clinical <- clinical[!is.na(clinical$response), ,drop=FALSE]
    x <- read_matrix(file.path(workspace,files[i]))
    if(i==1) rownames(x) <- hgnc$symbol[match(rownames(x),hgnc$entrez_id)]
    xs[[i]] <- rank_features(x[,rownames(clinical),drop=FALSE])
    ys[[i]] <- clinical$response
}
train <- do.call(rbind,xs); y <- unlist(ys)
mu <- colMeans(train); sdv <- apply(train,2,sd)
train <- scale(train,center=mu,scale=sdv)
test <- scale(rank_features(xp),center=mu,scale=sdv)
frozen <- jsonlite::fromJSON(file.path(output,"frozen_models.json"),simplifyVector=FALSE)
models <- c("IFN4_rank","IFN4_APM6_rank")
transport <- read.csv(file.path(output,"PUCH_frozen_model_performance.csv"),
                      row.names=1,check.names=FALSE)
coef_error <- c(); brier_errors <- c(); external_auc_errors <- c()
for(i in seq_along(models)) {
    m <- cbind(1, train[,seq_len(i),drop=FALSE])
    objective <- function(b) sum(log1p(exp(m%*%b)) - y*(m%*%b)) + sum(b[-1]^2)/2
    grad <- function(b) as.vector(crossprod(m,plogis(m%*%b)-y)) + c(0,b[-1])
    hess <- function(b) {
        pr <- as.vector(plogis(m%*%b))
        crossprod(m,(pr*(1-pr))*m) + diag(c(0,rep(1,ncol(m)-1)))
    }
    estimate <- nlminb(rep(0,ncol(m)),objective,gradient=grad,hessian=hess,
                       control=list(rel.tol=1e-14,x.tol=1e-14))
    error <- max(abs(grad(estimate$par)))
    stopifnot(error < 1e-5)
    spec <- frozen$models[[models[i]]]
    coef_error[i] <- max(abs(estimate$par - c(spec$intercept,unlist(spec$slopes))))
    preds <- plogis(cbind(1,test[,seq_len(i),drop=FALSE]) %*% estimate$par)
    external_auc_errors[i] <- abs(rank_auc(preds,yext)-transport[models[i],"auc"])
    brier_errors[i] <- abs(mean((yext-preds)^2)-transport[models[i],"brier"])
}
stopifnot(max(diff_auc,diff_rho,external_auc_errors)<1e-12, max(brier_errors)<1e-7,
          max(coef_error)<1e-5)
report <- list(status="passed", tools=R.version.string,
               patient_populations=c(training=length(y),external=length(yext)),
               standalone_rows_checked=nrow(stand), standalone_auc_max_abs_error=max(diff_auc),
               spearman_max_abs_error=max(diff_rho), frozen_coefficient_max_abs_error=max(coef_error),
               external_frozen_auc_max_abs_error=max(external_auc_errors),
               external_frozen_brier_max_abs_error=max(brier_errors),
               scope="Independent raw-to-point-estimate replay (base R ranks/U-statistic and nlminb ridge optimization). Main runner separately checks previously archived two-score intervals; this R replay does not verify the new interval endpoints.")
jsonlite::write_json(report,file.path(output,"independent_R_verification.json"),
                     pretty=TRUE,auto_unbox=TRUE,digits=15)
print(report)
