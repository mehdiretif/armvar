#' Welch's ANOVA mean comparison test for a given cell
#'
#' @param cell The name of the cell of interest.
#' @param geneMetadata The annotated normalized and scaled expression matrix.
#'
#' @returns The metrics of the Welch's ANOVA comparison test.
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' data <- armvar_processing(data = gao_2021_breast_cancer_data)
#' geneMetadata <- generate_geneMetadata(data)
#' welch_anova_res <- run_welch_anova(cell = "AAACCTGCACCTTGTC_b" , geneMetadata)
run_welch_anova <- function(cell, geneMetadata){
  res <- stats::oneway.test(stats::as.formula(paste0("`", cell, "` ~ chr")), data = geneMetadata, var.equal = FALSE)
  pval <- res$p.value
  f <- res$statistic
  c(cell = cell, pval = pval, fstat = f)
}

#' Testing the presence of mean differences between chromosome arms for every cell
#'
#' @param data The Seurat object.
#' @param geneMetadata The annotated normalized and scaled expression matrix.
#' @param n.cores The number of cores.
#'
#' @returns The Seurat object with the metrics of the Welch's ANOVA mean comparison test for all the tested cells.
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' data <- armvar_processing(data = gao_2021_breast_cancer_data)
#' geneMetadata <- generate_geneMetadata(data)
#' stat_test(data, geneMetadata, n.cores = 1)
#'
stat_test <- function(data, geneMetadata, n.cores=3){
  print("----- Welch's ANOVA Mean comparison test -----")
  cells <- colnames(geneMetadata[, 4:ncol(geneMetadata)])
  res <- mcprogress::pmclapply(cells, run_welch_anova, geneMetadata = geneMetadata, mc.cores = n.cores)
  stat_res <- as.data.frame(do.call(rbind, res))
  stat_res$fdr <- stats::p.adjust(stat_res$pval, method="bonferroni")
  stat_res[["welch_anova"]] <- ifelse(stat_res$fdr <= 0.05, 'cancer','normal')
  row.names(stat_res) <- stat_res$cell
  data <- Seurat::AddMetaData(data, stat_res)
  return(data)
}
