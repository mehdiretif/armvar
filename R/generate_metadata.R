#' Determine the position of genes in chromosome arms.
#'
#' @param data The Seurat object containing the normalized and scaled expression matrix in the 'SCT' assay.
#'
#' @returns The annotated normalized and scaled matrix.
#' @export
#'
#' @examples
#' data("gao_2021_breast_cancer_data", package = "armvar")
#' data <- armvar_processing(data = gao_2021_breast_cancer_data)
#' geneMetadata <- generate_geneMetadata(data)
generate_geneMetadata <- function(data, organism='human'){
  print("----- Gene annotation -----")

  if(organism=='human'){
    data("geneMetadata", package = "armvar")
  }
  geneMetadata <- geneMetadata[which(
    geneMetadata$gene_biotype %in% c("protein_coding","lncRNA") &
      geneMetadata$chromosome_name %in% c(1:22,"X") &
      geneMetadata$hgnc_symbol !=""),]

  geneMetadata <- unique(geneMetadata[,c("hgnc_symbol", "chromosome_name", "chr_arm", 'start_position')])

  scaled_data <- data[['SCT']]$scale.data

  #geneMetadata filtering based on the genes in data
  geneMetadata <- geneMetadata[geneMetadata$hgnc_symbol %in% row.names(scaled_data),]
  geneMetadata$chr <- paste0(geneMetadata$chromosome_name, geneMetadata$chr_arm)

  #gene filtering based on genes in metadata
  scaled_data <- scaled_data[row.names(scaled_data) %in% geneMetadata$hgnc_symbol,]

  geneMetadata <- merge(x=geneMetadata[,c('hgnc_symbol', 'chr', 'start_position')], y=scaled_data, by.x='hgnc_symbol', by.y='row.names')

  chr_to_remove <- c(names(table(geneMetadata$chr)[table(geneMetadata$chr) < 50]))
  if(length(chr_to_remove) > 0){
    geneMetadata <- geneMetadata[!(geneMetadata$chr %in% chr_to_remove),]
  }
  return(geneMetadata)
}
