#' Simulation of cancer tissues
#'
#' @param data A Seurat object containing cancer data.
#' @param sample_nb_to_simulate The number of cancer sample to simulate.
#' @param cancer_cell_nb_to_simulate The number of cancer cells to simulate in each sample.
#' @param total_cell_in_sample The total number of cells (cancer cells, healthy counterpart of the malignant cells and microenvironment cells) constituting one sample.
#' @param healthy_cells_to_transform A vector containing the names of the healthy cells to transform into cancer cells by inducing CNAs. We recommend to use the healthy counterparts of the malignant cells of interest.
#' @param tumor_microenvrnt_cells A vector containing the names of the microenvironment cells to use for the simulation (the cells will be added without modification).
#' @param sample_fc_list The list of the fold changes determined in the tumor samples.
#' @param sample_col The name of the column specifying samples in the metadata (data@meta.data).
#' @param seed A seed used to draw the cells and the fold changes that will be used for the simulation.
#'
#' @returns A Seurat object of the simulated data.
#' @export
#'
#' @examples
#' data("peng_2019_pdac_data")
#' cells_to_keep <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$orig.ident != 'T12',])
#' peng_2019_pdac_data <- subset(peng_2019_pdac_data, cells = cells_to_keep)
#' healthy_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 1' & peng_2019_pdac_data@meta.data$tissue == 'h',])
#' cancer_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 2',])
#' foldchange_computation_res <- foldchange_computation(peng_2019_pdac_data, healthy_cells, cancer_cells, sample_col='orig.ident', topNGenes=15000)
#' tumor_microenvrnt_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$tissue == "pt" & peng_2019_pdac_data@meta.data$published_annotation != 'Ductal cell type 2',])
#' simulated_seurat_obj <- cancer_tissue_simulation(data = peng_2019_pdac_data, sample_nb_to_simulate=23, cancer_cell_nb_to_simulate=320, total_cell_in_sample=400, healthy_cells_to_transform=healthy_cells, tumor_microenvrnt_cells=tumor_microenvrnt_cells, sample_fc_list= foldchange_computation_res$sample_fc_list,  sample_col='orig.ident', seed=1)
#'
cancer_tissue_simulation <- function(data, sample_nb_to_simulate, cancer_cell_nb_to_simulate, total_cell_in_sample, healthy_cells_to_transform, tumor_microenvrnt_cells, sample_fc_list, sample_col, seed=1){

  # check of the parameters
  micro_envrt_cell_nb <- total_cell_in_sample - cancer_cell_nb_to_simulate
  total_micro_envrt_cell_nb <- micro_envrt_cell_nb*sample_nb_to_simulate

  if(cancer_cell_nb_to_simulate > total_cell_in_sample){
    stop(paste0('ERROR: more cancer cells within one sample (', cancer_cell_nb_to_simulate, ') than the total number of cell in the given cells (', total_cell_in_sample,')'))
  }

  if(total_micro_envrt_cell_nb > length(tumor_microenvrnt_cells)){
    stop(paste0('ERROR: not enough tumor microenvrnt cells (expect at least ', total_micro_envrt_cell_nb, ' microenvironment cells'))
  }
  if(sample_nb_to_simulate <= length(sample_fc_list)){
    set.seed(seed)
    selected_tumor_sample_vector <- sample(names(sample_fc_list), size = sample_nb_to_simulate, replace = FALSE)
  }else{
    stop(paste0('ERROR: more copy number patterns than there are tumor samples (please specify a number of tumor sample > ', length(sample_fc_list)))
  }

  print(paste0(sample_nb_to_simulate, ' tumor samples containing ', cancer_cell_nb_to_simulate, ' tumor cell(s) each will be generated.'))
  count_data <- data[['RNA']]$counts

  simulated_count_list <- lapply(1:sample_nb_to_simulate, function(s){

    if(micro_envrt_cell_nb > 0){
      # select microenvrt cells
      set.seed(seed)
      collected_microenvrt_cells <- sample(tumor_microenvrnt_cells, micro_envrt_cell_nb, replace=FALSE)
      tumor_microenvrnt_cells <<- tumor_microenvrnt_cells[!tumor_microenvrnt_cells %in% collected_microenvrt_cells]

      simulated_counts <- as.data.frame(count_data[, collected_microenvrt_cells])
      colnames(simulated_counts) <- paste0('S', s, '_', colnames(simulated_counts))
    }

    # if there are cancer cells to simulate
    if(cancer_cell_nb_to_simulate > 0){

      # select the healthy cells which will be transformed into malignant cells
      healthy_metadata <- data@meta.data[healthy_cells_to_transform,]
      collected_ctrl_cells <- collect_healthy_to_transform(healthy_metadata, cancer_cell_nb_to_simulate, sample_col, seed=seed+s) # seed + s to vary the tissue or origin

      # select the cancer cells to determine the foldchanges that will be used for the transformation
      selected_tumor_sample <- selected_tumor_sample_vector[s]
      cancer_simulated_cell_name <- sapply(1:cancer_cell_nb_to_simulate, function(i){paste0('SIMULATED_S', s, '_', names(collected_ctrl_cells[i]), '_', selected_tumor_sample, '_', i)})
      foldchange <- sample_fc_list[[selected_tumor_sample]]

      # multiplying the counts by the foldchanges (VERIFIER QUE CA FONCTIONNE AVEC LES RAW COUNTS)
      cancer_simulated_counts <- count_data[, collected_ctrl_cells]
      cancer_simulated_counts[names(foldchange),] <- cancer_simulated_counts[names(foldchange),] * foldchange
      colnames(cancer_simulated_counts) <- cancer_simulated_cell_name

      if(micro_envrt_cell_nb > 0){
        simulated_counts <- cbind(cancer_simulated_counts, simulated_counts)
      }else{
        simulated_counts <- cancer_simulated_counts
      }
    }
    simulated_counts
  })
  simulated_count_mtx <- do.call(cbind, simulated_count_list)
  seurat_obj <- seurat_obj_construction(data, simulated_count_mtx)
  return(seurat_obj)
}

#' Collect similar healthy cells that will be transformed into malignant cells.
#'
#' @param healthy_metadata The Seurat metadata of the healthy cells.
#' @param cancer_cell_nb The number of cancer cells to simulate in each sample.
#' @param sample_col The name of the column specifying samples in the metadata (data@meta.data).
#' @param seed A seed used to draw the healthy cells that will be used for the simulation.
#'
#' @returns A vector of similar healthy cells (from the same sample) that will be used for the simulation.
#' @export
#'
#' @examples
#' data("peng_2019_pdac_data")
#' healthy_cells <- row.names(peng_2019_pdac_data@meta.data[peng_2019_pdac_data@meta.data$published_annotation == 'Ductal cell type 1' & peng_2019_pdac_data@meta.data$tissue == 'h',])
#' healthy_metadata <- peng_2019_pdac_data@meta.data[healthy_cells_to_transform,]
#' collected_healthy_cells <- collect_healthy_to_transform(healthy_metadata, cancer_cell_nb=100, sample_col='orig.ident', seed=1)
collect_healthy_to_transform <- function(healthy_metadata, cancer_cell_nb, sample_col, seed=1){

  collected_healthy_cells <- c()
  ctrl_samples <- unique(healthy_metadata[,sample_col])

  while(length(collected_healthy_cells) < cancer_cell_nb){
    remaining_cells_to_select <- cancer_cell_nb - length(collected_healthy_cells)

    set.seed(seed)
    selected_ctrl_sample <- sample(ctrl_samples, size = 1)
    ctrl_cells <- row.names(healthy_metadata[healthy_metadata[,sample_col] == selected_ctrl_sample,])
    ctrl_samples <- ctrl_samples[!ctrl_samples %in% selected_ctrl_sample] # remove the selected sample from the list of remaining ctrl samples (for the next random drawing)

    if(length(ctrl_cells) < remaining_cells_to_select){ # if less ctrl cells than required have been collected
      selected_cells <- ctrl_cells
      selected_cells <- setNames(selected_cells, rep(selected_ctrl_sample, length(ctrl_cells)))
    }
    if(length(ctrl_cells) > remaining_cells_to_select){ # if more ctrl cells than required have been collected
      set.seed(seed)
      selected_cells <- sample(ctrl_cells, remaining_cells_to_select, replace = FALSE)
      selected_cells <- setNames(selected_cells, rep(selected_ctrl_sample, remaining_cells_to_select))
    }
    collected_healthy_cells <- c(collected_healthy_cells, selected_cells)
  }
  return(collected_healthy_cells)
}

#' Construction and annotation of the Seurat object from the simulated count matrix.
#'
#' @param data The initial Seurat object from which the simulated count matrix was generated.
#' @param simulated_count_mtx The simulated count matrix.
#'
#' @returns The annotated Seurat object for the simulated data.
#' @export
#'
seurat_obj_construction <- function(data, simulated_count_mtx){

  simulated_seurat_obj <- CreateSeuratObject(as.matrix(simulated_count_mtx), assay = 'RNA')
  simulated_seurat_obj@meta.data$published_annotation <- ifelse(simulated_seurat_obj@meta.data$orig.ident == 'SIMULATED', 'cancer (simulated)', data@meta.data[sub("^[^_]*_", "", row.names(simulated_seurat_obj@meta.data)), 'published_annotation'])

  simulated_seurat_obj@meta.data$tumor_sample_origin <- ifelse(simulated_seurat_obj@meta.data$orig.ident == "SIMULATED", sapply(strsplit(rownames(simulated_seurat_obj@meta.data), "_"), `[`, 4), "healthy cells")

  simulated_seurat_obj@meta.data$group_origin <- ifelse(simulated_seurat_obj@meta.data$orig.ident == "SIMULATED", sapply(strsplit(rownames(simulated_seurat_obj@meta.data), "_"), `[`, 3), "healthy cells")

  simulated_seurat_obj@meta.data$sample <-  ifelse(simulated_seurat_obj@meta.data$orig.ident == "SIMULATED",  sapply(strsplit(rownames(simulated_seurat_obj@meta.data), "_"), `[`, 2), sapply(strsplit(rownames(simulated_seurat_obj@meta.data), "_"), `[`, 1))

  return(simulated_seurat_obj)
}

