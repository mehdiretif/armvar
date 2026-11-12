
#' Smoothing of the Welch's ANOVA classification using nearest neighbors.
#'
#' @param data The Seurat object containing the scaled matrix.
#' @param k The number of nearest neighbors.
#' @param n_pcs The number of principal components.
#' @param n.cores The number of cores.
#'
#' @returns The Seurat object with the ArmVar final classification and smoothed F-stat (in the meta.data).
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' data <- armvar_processing(data = gao_2021_breast_cancer_data)
#' geneMetadata <- generate_geneMetadata(data)
#' data <- stat_test(data, geneMetadata, n.cores = 1)
#' armvar_res <- armvar_smoothing(data, k = c(200, 300), n_pcs = 50, n.cores = 1)
#'
armvar_smoothing <- function(data, k, n_pcs, n.cores){

  print('----- knn smoothing -----')
  data <- Seurat::RunPCA(data, verbose = FALSE, assay = 'SCT', seed.use = 42)
  cell_emb <- data@reductions$pca@cell.embeddings

  combi <- expand.grid(k, n_pcs)
  armvar_list <- lapply(1:nrow(combi), function(c){
    k <- combi[c, 1]
    n_pcs <- combi[c, 2]
    print(paste0('K=', k, ' - PCs=', n_pcs))

    data <- Seurat::FindNeighbors(data, reduction = "pca", dims = 1:n_pcs, verbose = FALSE, k.param = k, assay = 'SCT')
    knn_data <- data@graphs$SCT_nn

    # knn smoothing
    knn_res <- do.call(rbind, mcprogress::pmclapply(row.names(knn_data), mc.cores = n.cores, function(cell){

      neighbors <- names(knn_data[cell, knn_data[cell, ] == 1]) # considering the cell itself as its nearest neighbor (distance = 0)
      #neighbors <- neighbors[neighbors != cell] #remove the cell from its neighbors
      dist <-  abs(colSums(apply(cell_emb[neighbors,], 1, function(x){x - cell_emb[cell,]}))) # euclidean distance computation
      coef <- 1/(log10(dist + 1)+1)

      neighbors_annotation <- data@meta.data[neighbors, 'welch_anova', drop=FALSE]
      cancer_neighbors <- row.names(neighbors_annotation[neighbors_annotation$welch_anova == 'cancer', ,drop=FALSE])
      normal_neighbors <- row.names(neighbors_annotation[neighbors_annotation$welch_anova == 'normal', ,drop=FALSE])

      #basic_knn_correction <- ifelse(length(cancer_neighbors) > length(normal_neighbors), 'cancer', ifelse(length(cancer_neighbors) < length(normal_neighbors), 'normal', stat_res[cell, paste0('annotation_', test_name)]))
      smoothened_class <- ifelse(sum(coef[cancer_neighbors]) > sum(coef[normal_neighbors]), 'cancer', 'normal')

      smoothened_fstat <- sum(as.numeric(data@meta.data[neighbors, 'fstat.F']) * coef)/sum(coef)
      stats::setNames(data.frame(cell, smoothened_class, smoothened_fstat), c('cell', paste0('ArmVar_k', k, '_pc', n_pcs), paste0('smoothed_fstat_k', k ,'_pc', n_pcs)))
    }))
  })
  if(length(armvar_list)>1){
    armvar_res <- Reduce(function(x, y) merge(x, y, by = "cell"), armvar_list)
  }else{
    armvar_res <- armvar_list[[1]]
  }
  row.names(armvar_res) <- armvar_res$cell
  data <- Seurat::AddMetaData(data, armvar_res)

 return(data)
}
