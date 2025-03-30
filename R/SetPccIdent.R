#' @title Set the ident as the default in a plantcellchat object.
#'
#' @description This function sets the specified `input.ident` as the default
#'     label for the PlantCellChat object. It can be set either from the meta
#'         information, if not found, as a cell label when the number of
#'             unique cell labels matches the number of columns in the data matrix.
#'
#' @param pcc_obj A PlantCellChat object.
#' @param input.ident The identifier to be set as the default label.
#'
#' @return A plantcellchat object after setting the default ident.
#'
#' @export
#'
#' @examples
#' \donttest{
#' pcc_obj <- SetPccIdent(pcc_obj, input.ident = "cell_type")
#' }
SetPccIdent <- function(pcc_obj,input.ident) {
  if (is.character(input.ident) && length(input.ident) == 1) {
    if (input.ident %in% colnames(pcc_obj@meta)) {
      pcc_obj@idents <- factor(pcc_obj@meta[[input.ident]])
      message(input.ident," was successfully set as the default label.",'\n')
    } else {
      stop("This ident is not found in the meta data.",'\n')
    }
  }

  if (!is.character(input.ident)) {
    if (length(input.ident) == ncol(pcc_obj@data)) {
      pcc_obj@idents <- factor(input.ident)
    } else {
      message("The length of input.ident is not equal to the number of columns in the sparse matrix.",'\n')
    }
  }

  return(pcc_obj)
}
