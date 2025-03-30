#' @title Create a shortcut command to print PlantCellChat ident.
#'
#' @param pcc_obj A PlantCellChat object.
#'
#' @return The PlantCellChat defalut ident.
#' @export
#'
#' @examples PccIdents(pcc_obj)
PccIdents <- function(pcc_obj){
  if (is.null(pcc_obj@idents)) {
    stop("No information available in 'idents' slot, please run `SetPccIdent` first.")
  }
  idents <- pcc_obj@idents
  return(idents)
}
