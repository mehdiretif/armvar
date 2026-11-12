#' Execute ArmVar
#'
#' @param data A Seurat object.
#' @param k The number of nearest neighbor cells to take into account for the smoothing of the initial classification (Welch's ANOVA).
#' @param n_pcs The number of principal components (PCs) to use to identify the nearest neighbors.
#' @param n.cores The number of cores.
#' @param heatmap TRUE to generate an heatmap to visualize the ArmVar comparison test (else FALSE).
#' @param heatmap_var The name of the columns in the data@meta.data to visualize in the heatmap legend ('sample' and 'published_annotation' by default).
#'
#' @returns The initial data with the results of ArmVar in the meta.data (Welch's ANOVA and ArmVar classifications, and the smoothed F-stat).
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' armvar(data = gao_2021_breast_cancer_data,
#'        k = 300,
#'        n_pcs = 50,
#'        n.cores = 1,
#'        heatmap = TRUE,
#'        heatmap_var = c('sample', 'published_annotation'))
#'
armvar <- function(data, k, n_pcs, n.cores = 5, heatmap = TRUE, heatmap_var = c('sample', 'published_annotation')){

  data <- armvar_processing(data)
  geneMetadata <- generate_geneMetadata(data) # generate the chromosome location associated with the genes

  data <- stat_test(data, geneMetadata, n.cores) # statistical test (Welch's ANOVA)

  data <- armvar_smoothing(data, k, n_pcs, n.cores)

  armvar_object <- list(data = data, geneMetadata = geneMetadata)

  if(heatmap == TRUE){
    hm_object <- generate_heatmap(data, geneMetadata, heatmap_var)
    armvar_object[['heatmap_obj']] <- hm_object
  }

  return(armvar_object)
}
