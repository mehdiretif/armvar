#' Heatmap construction to visualize average differences between chromosome arms
#'
#' @param data The Seurat object with the ArmVar results in the meta.data.
#' @param geneMetadata The annotated normalized and scaled expression matrix.
#' @param heatmap_var The name of the columns in the meta.data that will be added to the heatmap legend ('sample' and 'published_annotation' by default).
#'
#' @returns A heatmap with the ArmVar results.
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' data <- armvar_processing(data = gao_2021_breast_cancer_data)
#' geneMetadata <- generate_geneMetadata(data)
#' data <- stat_test(data, geneMetadata, n.cores = 1)
#' data <- armvar_smoothing(data, k = c(200, 300), n_pcs = 50, n.cores = 5)
#' heatmap_obj <- generate_heatmap(data = data,
#'                  geneMetadata = geneMetadata,
#'                  heatmap_var = c('sample', 'published_annotation'))
#' ComplexHeatmap::draw(heatmap_obj$heatmap, annotation_legend_list = heatmap_obj$legend)
generate_heatmap <- function(data, geneMetadata, heatmap_var = c('sample', 'published_annotation')){

  print(paste0('----- heatmap generation -----'))

  chr_order <- paste0(rep(c(1:22, "X"), each = 2), rep(c("p", "q"), times = 23))
  chr_order <- intersect(chr_order, unique(geneMetadata$chr))

  scaled_data <- data[['SCT']]$scale.data

  chr_average_list <- stats::setNames(lapply(chr_order, function(chr){
    genes <- geneMetadata[geneMetadata$chr == chr, 'hgnc_symbol']
    colMeans(scaled_data[genes, , drop=FALSE])
  }), chr_order)
  hm_mtx <- do.call(rbind, chr_average_list)

  # cell ordering depending of the given var
  sorted_cells <- row.names(data@meta.data[do.call(order, data@meta.data[, heatmap_var, drop = FALSE]),])
  hm_mtx <- hm_mtx[, sorted_cells]

  chr_annotation = ComplexHeatmap::HeatmapAnnotation(
    chr = ComplexHeatmap::anno_block(
      gp = grid::gpar(
        fill = rep(
          c("gray26", "gray61"),
          length.out = nrow(hm_mtx)
        ),
        col = "black",
        lwd = 0.3
      )
    ),
    annotation_height = grid::unit(6, "mm")
  )

  armvar_col <- colnames(data@meta.data)[grepl('ArmVar', colnames(data@meta.data))]
  armvar_parameters <- gsub("ArmVar_", '', armvar_col)
  armvar_var <- unlist(lapply(armvar_parameters, function(p){colnames(data@meta.data)[grepl(p, colnames(data@meta.data))]}))

  fstat_col <- colnames(data@meta.data)[grepl('smoothed_fstat', colnames(data@meta.data))]
  legend_var <- c(heatmap_var, armvar_var)

  legend_value_list <- stats::setNames(lapply(legend_var, function(v){data@meta.data[sorted_cells, v]}), vname(legend_var))
  var_colors <- color_selection(data, heatmap_var)

  # legend bars (on the left of the heatmap)
  armvar_colors <- stats::setNames(lapply(armvar_var, function(col){
    if(grepl('ArmVar', col)){
      c("normal" = "#00BFC4", "cancer" = "#F8766D")
    }else{
      circlize::colorRamp2(c(0, ceiling(max(data@meta.data[, fstat_col]))), c("white", "blue"))
    }
  }), vname(armvar_var))

  colors <- c(armvar_colors, var_colors)

  cell_annotation = do.call(ComplexHeatmap::rowAnnotation,
                            c(legend_value_list,
                              list(
                                show_legend = FALSE,
                                col = colors
                              )
                            )
  )

  # legend
  armvar_legend = ComplexHeatmap::Legend(title = "ArmVar", at = c("normal", "cancer"), legend_gp = grid::gpar(fill = c("#00BFC4", "#F8766D")))
  fstat_legend = ComplexHeatmap::Legend(title = "Smoothed F-stat", col_fun = circlize::colorRamp2(c(0, ceiling(max(data@meta.data[, fstat_col]))), c("white", "blue")))
  heatmap_var_legend <- lapply(names(var_colors), function(col){
    ComplexHeatmap::Legend(title = col, at = names(var_colors[[col]]), legend_gp = grid::gpar(fill = c(var_colors[[col]])))
  })

  legend <- do.call(ComplexHeatmap::packLegend,
                    c(list(armvar_legend, fstat_legend),
                      unlist(heatmap_var_legend),
                      list(direction = "vertical", max_height = grid::unit(8, "cm")))
  )

  hm_mtx <- sign(hm_mtx) * log1p(abs(hm_mtx))
  scale_val <- round(max(max(hm_mtx), abs(min(hm_mtx))),2)

  heatmap_legend = list(title =  "log1p(average(scaled data))", at = c(-scale_val,0,scale_val), labels = c(paste0(-scale_val, " (potential loss)"), "0", paste0(scale_val, " (potential gain)")))

  hm <- ComplexHeatmap::Heatmap(
    matrix = t(hm_mtx),
    show_heatmap_legend = TRUE,
    heatmap_legend_param = heatmap_legend,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    top_annotation = chr_annotation,
    left_annotation = cell_annotation,
    use_raster = FALSE,
    column_split = row.names(hm_mtx),
    column_gap = grid::unit(0, "mm"),
    column_title = 'chromosome arms',
    column_title_gp = grid::gpar(fontsize = 12),
  )
  hm_object <- list(heatmap = hm, legend = legend)

  return(hm_object)
}

#' Select a ggtheme palette among 'Winter', 'Tableau_10', 'Tableau_20', 'Superfishel_Stone', 'Miller_Stone', 'manyeys', 'Hue_Circle', 'Classic_10_Medium', 'gdoc' to visualize the specified column(s).
#'
#' @param data The Seurat object
#' @param heatmap_var The column(s) to visualize in the heatmap legend.
#'
#' @returns A list with color vectors for the specified columns.
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' color_res <- color_selection(data, 'published_annotation')
color_selection <- function(data, heatmap_var){
  var_colors <- list()
  palette_table <- paletteer::palettes_d_names
  palettes <- list('ggthemes' = c('Winter', 'Tableau_10', 'Tableau_20', 'Superfishel_Stone', 'Miller_Stone', 'manyeys', 'Hue_Circle', 'Classic_10_Medium', 'gdoc'),
                   'colorBlindness' = c('SteppedSequential5Steps'),
                   'khroma' = c('soil', 'discreterainbow'),
                   'Polychrome' = c('light', 'sky', 'palette36'))
  palette_table <- palette_table[palette_table$palette %in% unlist(palettes) & palette_table$package %in% names(palettes),]
  for(v in heatmap_var){
    values <- sort(unique(data@meta.data[,v]))
    size <- length(values)
    difference <- stats::setNames(palette_table$length - size, palette_table$palette)
    set.seed(1)
    selected_palette <- sample(names(difference[difference == min(difference[difference>0])]), 1)
    palette_table <- palette_table[palette_table$palette != selected_palette,]

    pac <- unlist(sapply(names(palettes), function(pac){if(selected_palette %in% palettes[[pac]]){pac}}))
    colors <- paletteer::paletteer_d(paste0(pac,"::", selected_palette))
    var_colors[[vname(v)]] <- stats::setNames(as.character(colors[1:size]), values)
  }
  return(var_colors)
}

#' Rename the variable for the heatmap legend.
#'
#' @param name The raw name.
#'
#' @returns The modified name.
#' @export
#'
#' @examples vname('ArmVar_k300_pc50')
vname <- function(name){
  name <- gsub('_', ' ', name)
  name <- sub(" (k\\d+) (pc\\d+)$", " (\\1, \\2)", name)
  name <- paste0(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)))
  return(name)
}
