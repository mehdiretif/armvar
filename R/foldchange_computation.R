#' Computation of fold changes between cancer and healthy cells for genes detected in chromosome arms with CNA(s) according to fastCNV (supervised)
#'
#' @param data A Seurat object containing cancer data
#' @param healthy_cells A vector of cells defined healthy by the user.
#' @param cancer_cells  A vector of tumor cells.
#' @param sample_col The name of the column specifying samples in the metadata (data@meta.data).
#' @param topNGenes The number of top expressed genes used by fastCNV.
#'
#' @returns
#' A list containing fold changes for the genes located in the chromosome arms with CNAs (according to fastCNV) computed for every cell with such alteration(s) (cell_fc_list).
#' The patterns of CNAs shared by the cells (patterns).
#' A list containing the accumulated fold changes for every sample. For genes associated with multiple fold changes (from multiple cells), the corresponding final fold change is determined using the median of the fold changes. (sample_fc_list).
#' The metadata with the cnv fraction scores computed by fastCNV (fastCNV_metadata).
#' @export
#'
#' @examples
#' data("peng_2019_pdac_data")
#' cells_to_keep <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$orig.ident != 'T12',])
#' peng_2019_pdac_data <- subset(peng_2019_pdac_data, cells = cells_to_keep)
#' healthy_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 1' & peng_2019_pdac_data@meta.data$tissue == 'h',])
#' cancer_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 2',])
#' foldchange_computation_res <- foldchange_computation(peng_2019_pdac_data, healthy_cells, cancer_cells, sample_col='orig.ident', topNGenes=15000)

foldchange_computation <- function(data, healthy_cells, cancer_cells, sample_col, topNGenes=15000){
  # filter the data by the cells of interest
  filtered_data <- subset(data, cells = c(healthy_cells, cancer_cells))
  filtered_data@meta.data$ground_truth <- ifelse(row.names(filtered_data@meta.data) %in% healthy_cells, 'reference', 'cancer')

  fastcnv_computation_res <- chr_scale_fastcnv(filtered_data, sample_col, topNGenes)

  CFC_res <- cell_foldchange_computation(filtered_data, healthy_cells, cancer_cells, sample_col, fastcnv_computation_res$genomicScores, fastcnv_computation_res$genomicWindows)

  sample_fc_list <- sample_FCmedian_calculation(filtered_data, CFC_res$foldchanges)

  return(list('cell_fc_list' = CFC_res$foldchanges, 'pattern_list' = CFC_res$patterns, 'sample_fc_list' = sample_fc_list, 'fastCNV_metadata' = fastcnv_computation_res$metadata))
}
#'
#' Use fastCNV to identify the chromosome arms containing CNA(s)
#'
#' @param data A Seurat object containing cancer data.
#' @param sample_col The name of the column specifying samples in the metadata (data@meta.data).
#' @param topNGenes The number of top expressed genes used by fastCNV.
#' @param windowStep The step of the window scanning the genome used by fastCNV to detect CNAs.
#' @param windowSize The size of the window scanning the genome used by fastCNV to detect CNAs.
#'
#' @returns
#' A list containing data specifying the genes located within the fastCNV windows (genomicWindows),
#' the fastCNV genomic scores (genomicScores),
#' and metadata with the cnv fraction scores computed by fastCNV (fastCNV_metadata).
#' @export
#'
#' @examples
#' data("peng_2019_pdac_data")
#' cells_to_keep <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$orig.ident != 'T12',])
#' filtered_data <- subset(peng_2019_pdac_data, cells = cells_to_keep)
#' fastcnv_computation_res <- chr_scale_fastcnv(filtered_data, sample_col='orig.ident', topNGenes=15000)
chr_scale_fastcnv <- function(data, sample_col, topNGenes, windowStep = 50, windowSize = 10000){

  print('#### fastCNV to detect altered regions ####')
  seurat_obj_list <- Seurat::SplitObject(data, split.by = sample_col)

  fastcnv_res <- fastCNV::fastCNV(seuratObj=seurat_obj_list, sampleName=names(seurat_obj_list), referenceVar = 'ground_truth', referenceLabel = rep('reference', length(seurat_obj_list)), scaleOnReferenceLabel = TRUE, prepareCounts = "default", getCNVPerChromosomeArm = FALSE, getCNVClusters = FALSE, doPlot = FALSE, topNGenes = topNGenes, windowSize=windowSize, windowStep=windowStep, saveGenomicWindows = TRUE)

  # matrix with average expression for duplicated or deleted regions (fastCNV)
  genomicScores <- as.matrix(do.call(cbind, lapply(fastcnv_res, function(x) x@assays$genomicScores$counts)))

  metadata <- do.call(rbind, lapply(fastcnv_res, function(x) x@meta.data))

  load(paste0("genomicWindows_size", windowSize, "_step", windowStep, ".RData"))
  return(list('genomicWindows' = genomicWindows, 'genomicScores' = genomicScores, 'metadata' = metadata))
}
#'
#' Computation of fold changes for the genes located in the chromosome arms with CNAs according to fastCNV
#'
#' @param data A Seurat object containing cancer data.
#' @param healthy_cells A vector of cells defined healthy by the user.
#' @param cancer_cells A vector of tumor cells.
#' @param sample_col The name of the column specifying samples in the metadata (data@meta.data).
#' @param genomicScores A matrix (windows x cells) containing the fastCNV genomic scores.
#' @param genomicWindows The list specifying the genes located within the fastCNV windows.
#'
#' @returns
#' A list containing fold changes for the genes located in the chromosome arms with CNAs (according to fastCNV) computed for every cell with such alteration(s) (foldchanges),
#' and the patterns of CNAs shared by the cells (patterns).
#' @export
#'
#' @examples
#' data("peng_2019_pdac_data")
#' cells_to_keep <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$orig.ident != 'T12',])
#' peng_2019_pdac_data <- subset(peng_2019_pdac_data, cells = cells_to_keep)
#' filtered_data <- subset(peng_2019_pdac_data, cells = c(healthy_cells, cancer_cells))
#' filtered_data@meta.data$ground_truth <- ifelse(row.names(filtered_data@meta.data) %in% healthy_cells, 'reference', 'cancer')
#' fastcnv_computation_res <- chr_scale_fastcnv(filtered_data, sample_col='orig.ident', topNGenes=15000)
#' healthy_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 1' & peng_2019_pdac_data@meta.data$tissue == 'h',])
#' cancer_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 2',])
#' CFC_res <- cell_foldchange_computation(filtered_data, healthy_cells, cancer_cells, sample_col='orig.ident', fastcnv_computation_res$genomicScores, fastcnv_computation_res$genomicWindows)
#'
cell_foldchange_computation <- function(data, healthy_cells, cancer_cells, sample_col, genomicScores, genomicWindows){

  # normalization without log transformation
  data <- Seurat::NormalizeData(data, assay = 'RNA', normalization.method = 'RC')
  norm_count_matrix <- as.matrix(data[['RNA']]$data)

  print('#### foldchange computation for every cells ####')
  complete_fc_list <- list()
  cell_list <- list()

  tumor_samples <- unique(data@meta.data[intersect(colnames(data), cancer_cells), sample_col])
  for(sample in tumor_samples){
    print(paste0('##### ', sample,' #####'))

    complete_fc_list[[sample]] <- list()
    cell_list[[sample]] <- list()

    tumor_cells <- row.names(data@meta.data[data@meta.data[, sample_col] == sample,])
    # filter the fastcnv data
    genomicScores_i <- genomicScores[, tumor_cells, drop = FALSE]
    # keep cells from tumor sample with cnvs
    wthcnv_genomicScores_i <- genomicScores_i[, colSums(genomicScores_i) != 0, drop = FALSE]
    wthcnv_genomicScores_i <- ifelse(wthcnv_genomicScores_i > 0, paste0("+_",row.names(wthcnv_genomicScores_i)), ifelse(wthcnv_genomicScores_i < 0, paste0("-_",row.names(wthcnv_genomicScores_i)), wthcnv_genomicScores_i))

    # collect all the CNVs patterns
    patterns <- apply(wthcnv_genomicScores_i, 2, paste0, collapse = "_")
    pattern_list <- split(colnames(wthcnv_genomicScores_i), patterns)

    complete_fc_list[[sample]] <- lapply(names(pattern_list), function(p) {
      splited_p <- strsplit(p, '_')[[1]]
      window_list <- splited_p[!(splited_p %in% c('0', '+', '-'))]

      cnv_genes <- unique(unlist(genomicWindows[window_list]))
      cnv_cell_list <- unlist(pattern_list[p])

      cnv_mean <- rowMeans(norm_count_matrix[cnv_genes, cnv_cell_list, drop = FALSE])
      ref_mean <- rowMeans(norm_count_matrix[cnv_genes, healthy_cells, drop = FALSE])
      fc_list <- (cnv_mean + 1e-9)/ (ref_mean + 1e-9)

      list(
        median = stats::median(fc_list),
        log2 = log2(fc_list),
        values = fc_list
      )
    })
    cell_list[[sample]] <- pattern_list
  }
return(list('foldchanges' = complete_fc_list, 'patterns' = cell_list))
}

#' Accumulation of cell fold changes for the tumor samples
#'
#' @param data A Seurat object containing cancer data.
#' @param complete_fc_list A list containing the fold changes computed for every cell.
#'
#' @returns A list containing the accumulated fold changes for every sample. For genes associated with multiple fold changes (from multiple cells), the corresponding final fold change is determined using the median of the fold changes.
#' @export
#'
#' @examples
#' data("peng_2019_pdac_data")
#' cells_to_keep <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$orig.ident != 'T12',])
#' peng_2019_pdac_data <- subset(peng_2019_pdac_data, cells = cells_to_keep)
#' filtered_data <- subset(peng_2019_pdac_data, cells = c(healthy_cells, cancer_cells))
#' filtered_data@meta.data$ground_truth <- ifelse(row.names(filtered_data@meta.data) %in% healthy_cells, 'reference', 'cancer')
#' fastcnv_computation_res <- chr_scale_fastcnv(filtered_data, sample_col='orig.ident', topNGenes=15000)
#' healthy_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 1' & peng_2019_pdac_data@meta.data$tissue == 'h',])
#' cancer_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 2',])
#' CFC_res <- cell_foldchange_computation(filtered_data, healthy_cells, cancer_cells, sample_col='orig.ident', fastcnv_computation_res$genomicScores, fastcnv_computation_res$genomicWindows)
#' sample_fc_list <- sample_FCmedian_calculation(filtered_data, CFC_res$foldchanges)
#'
sample_FCmedian_calculation <- function(data, complete_fc_list){

  print('foldchange accumulation')
  median_fc_list <- list()
  for(sample in names(complete_fc_list)){
    fc_df <- data.frame(row.names = row.names(data))
    for(cell in 1:length(complete_fc_list[[sample]])){
      fc_df[names(complete_fc_list[[sample]][[cell]]$values), paste('cell_', cell)] <- complete_fc_list[[sample]][[cell]]$values
    }
    # median of the gene foldchanges of the cells belonging to the same sample (without accounting the NA)
    median_fc_list[[sample]] <- setNames(matrixStats::rowMedians(as.matrix(fc_df), na.rm = TRUE), row.names(fc_df))
    median_fc_list[[sample]][is.na(median_fc_list[[sample]])] <- 1
  }
  return(median_fc_list)
}


