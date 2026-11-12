#' Processing
#'
#' @param data A Seurat object.
#'
#' @returns The Seurat object with the normalized and scaled expression matrix in the SCT assay.
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' armvar_processing(data = gao_2021_breast_cancer_data)
#'
armvar_processing <- function(data){
  print('----- Preprocessing step -----')
  genes_to_keep <- row.names(data[['RNA']]$counts)[Matrix::rowSums(data[['RNA']]$counts > 0)/ncol(data[['RNA']]$counts) >= 0.05]
  HLA_genes <- row.names(data[['RNA']]$counts)[grep('HLA-', row.names(data[['RNA']]$counts))]
  cell_cycle_genes <- unlist(Seurat::cc.genes)
  genes_to_keep <- genes_to_keep[!(genes_to_keep %in% c(HLA_genes, cell_cycle_genes))]
  print(paste0(nrow(data[['RNA']]$counts) - length(genes_to_keep), ' genes removed'))
  data <- subset(data, features = genes_to_keep)

  print('----- Normalization and Scaling (SCTransform) -----')
  data <- Seurat::SCTransform(data, vars.to.regress = 'nFeature_RNA', verbose = FALSE,  do.scale = FALSE, do.center = TRUE, return.only.var.genes = FALSE, seed.use = 1448145)

  return(data)
}
