#' @title Extract subset sparse matrix of ligands and receptors gene which mediating signaling.
#'
#' @param pcc_obj A PlantCellChat object.
#'
#' @return A sparse matrix recording ligands and receptors for signaling which saved in
#'     the `data.signaling` slot of the PlantCellChat object.
#' @export
#'
ExtractSignalingData <- function(pcc_obj) {
  data <- pcc_obj@data
  genes <- pcc_obj@database$genelist$Gene_ID

  if (!any(genes %in% rownames(data))) {
    stop("Ligands and receptors do not exist in the sparse matrix, please check whether the correct database is being used.")
  }

  pcc_obj@data.signaling <- data[rownames(data) %in% genes,]
  message(length(intersect(rownames(pcc_obj@data),genes)), " genes have been extracted. The new sparse matrix saved in the 'data.signaling' slot.")

  return(pcc_obj)
}
