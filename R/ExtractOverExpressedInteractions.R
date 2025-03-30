#' @title Extract over-expressed interactions
#'
#' @description
#' Extract over-expressed ligand-receptor interactions from PlantCellChatDB.
#'
#' @param pcc_obj A PlantCellChat object which contains the result of differential expression analysis(saved in the 'diffexp' slot).
#'
#' @return The over-expressed ligand-receptor interactions which saved in
#'     the `LR` slot of the PlantCellChat object.
#'
#' @export
#'
ExtractOverExpressedInteractions <- function(pcc_obj) {
  interactions <- pcc_obj@database$interaction
  genes <- pcc_obj@diffexp$features

  # Extracting signaling molecules and receptors interactions
  SRI <- subset(interactions,
                Interaction_type == "Signaling molecule" &
                  Receptor %in% genes)

  # Extracting ligands and receptors interactions
  LRI <- subset(interactions,
                Interaction_type == "Interacting protein" &
                  Receptor %in% genes &
                  Ligand %in% genes)

  # Merging all interactions
  DE_interactions <- rbind(SRI,LRI)
  rownames(DE_interactions) <- DE_interactions$Interaction_name
  pcc_obj@lrpairs$lr.signaling <- DE_interactions

  message(paste0(nrow(DE_interactions), " over-expressed ligand-receptor interactions have been extracted. The ligand-receptor pairs information saved in the `LR` slot."))
  return(pcc_obj)
}
