#' @title Adjust Communication Strength Based on Spatial Distance
#'
#' @description This function adjusts the communication strength between cells based on their spatial distance.
#'
#' @param prob A 3D array representing communication probabilities between cells.
#' @param coord A matrix of spatial coordinates for each cell. Rows are cell names, and columns are spatial dimensions.
#' @param neighbors A named list where each element contains the nearest neighbors of a cell, as determined by a k-nearest neighbors algorithm.
#' @param n.distance Numeric, the spatial distance constant for adjusting communication strength.
#'
#' @return A 3D array representing communication probabilities between cells, adjusted based on spatial distance.
#' @export

AdjustCommunStrength <- function(prob, coord, neighbors, n.distance) {
  for (cell1 in rownames(coord)) {
    current_coord <- coord[cell1, ]

    for (cell2 in rownames(coord)) {
      if (cell2 != cell1 && !(cell2 %in% neighbors[[cell1]])) {

        distance <- sqrt(sum((current_coord - coord[cell2, ])^2))
        factor <- n.distance / distance
        prob[cell1, cell2, ] <- prob[cell1, cell2, ] * factor
      }
    }
  }
  return(prob)
}

#' @title Calculate Spatial Communication Strength
#'
#' @description This function calculates the spatial communication strength between cells based on ligand-receptor interactions.
#'
#' @param pcc_obj A PlantCellChat object
#' @param Kh Numeric, the dissociation constant in the Hill equation, which modulates the sensitivity of signaling molecule binding.
#' @param n Numeric, the Hill coefficient, which describes the cooperativity of signaling molecule binding.
#' @param k.neighbors Integer, the number of nearest neighbors to consider in spatial calculations. Default is 5.
#' @param n.distance Numeric, the spatial distance constant for adjusting communication strength. Default is 10.
#'
#' @return A PlantCellChat object with updated spatial communication strength information stored in the `sp.net` slot.
#' @export

CalculateSpatialCommunStrength <- function(pcc_obj,Kh = 0.5, n = 1, coord, k.neighbors = 5, n.distance = 10) {
  if (is.null(pcc_obj@lrpairs$lr.signaling)) {
    stop("Over-expressed interaction were not extracted, please run the function 'ExtractOverExpressedInteractions' first.")
  } else {
    interaction <- pcc_obj@lrpairs$lr.signaling
    interaction <- interaction[interaction$Interaction_type == "Interacting protein",]
    if (is.null(pcc_obj@data.signaling)) {
      stop("The signaling data were not extracted, please run the function 'ExtractSignalingData' first.")
    } else {
      exp <- as.matrix(pcc_obj@data.signaling)
    }
  }

  # Ligand
  exp_L <- ExtractExpMatrix(pcc_obj,comps = "Ligand",interaction.type = "lri",na.replace = T,spatial = T)

  # Receptor
  exp_single_R <- ExtractExpMatrix(pcc_obj,comps = "Receptor",interaction.type = "lri",na.replace = T,spatial = T)

  exp_complex1 <- ExtractExpMatrix(pcc_obj,comps = "Complex1",interaction.type = "lri",spatial = T)
  exp_complex2 <- ExtractExpMatrix(pcc_obj,comps = "Complex2",interaction.type = "lri",spatial = T)
  sub1unit <- exp_single_R * correct(exp_complex1)
  sub2unit <- (exp_single_R * exp_complex1 * correct(exp_complex2))^(1/2)
  sub3unit <- (exp_single_R * exp_complex1 * exp_complex2)^(1/3)

  exp_coI <- ExtractExpMatrix(pcc_obj,comps = "CoAreceptor",interaction.type = "lri",zero.replace = T,spatial = T)
  exp_coA <- ExtractExpMatrix(pcc_obj,comps = "CoIreceptor",interaction.type = "lri",zero.replace = T,spatial = T)

  exp_R <- MergeExpMatrix(matrix_list = list(sub1unit,sub2unit,sub3unit)) * (1 + exp_coA) * (1/(1 + exp_coI))

  # Agonist
  exp_ag <- ExtractExpMatrix(pcc_obj,comps = "Agonist",interaction.type = "lri",zero.replace = T,spatial = T)
  p_ag <- 1 + (exp_ag^n/(Kh^n + exp_ag^n))

  # Antagonist
  exp_anta <- ExtractExpMatrix(pcc_obj,comps = "Antagonist",interaction.type = "lri",zero.replace = T,spatial = T)
  p_anta <- Kh^n / (Kh^n + exp_anta^n)

  # calculate communication strength
  nLR <- nrow(interaction)
  numcell <- ncol(exp)
  prob <- array(NA,dim = c(numcell,numcell,nLR))

  for (i in 1:nLR){
    prob[,,i] <- (exp_L[i,] %*% t(exp_R[i,]))^n/(Kh^n + (exp_L[i,] %*% t(exp_R[i,]))^n) * (p_ag[i,] %*% t(p_ag[i,])) * (p_anta[i,] %*% t(p_anta[i,]))
  }

  dimnames(prob)[[1]] <- colnames(exp)
  dimnames(prob)[[2]] <- colnames(exp)
  dimnames(prob)[[3]] <- rownames(interaction)

  coord <- pcc_obj@spatial$cell.embeddings

  knn <- FNN::get.knn(coord, k = k.neighbors)

  neighbors <- lapply(1:nrow(knn$nn.index), function(i) {
    rownames(coord)[knn$nn.index[i, ]]
  })

  names(neighbors) <- rownames(coord)

  adj_prob <- AdjustStrength(prob, coord, neighbors, n.distance = n.distance)

  df <- data.frame()
  for (i in 1:dim(adj_prob)[[3]]) {
    source_name <- rownames(adj_prob)[i]
    df0 <- as.data.frame(adj_prob[,,i]) %>%
      tibble::rownames_to_column(var = "Target") %>%
      tidyr::pivot_longer(cols = -Target, names_to = "Source", values_to = "Prob") %>%
      dplyr::filter(Prob != 0) %>%
      dplyr::mutate(Source = source_name,
                    Interaction_name = dimnames(adj_prob)[[3]][i]) %>%
      dplyr::left_join(interaction,by = "Interaction_name") %>%
      dplyr::select(Interaction_name, Signal, Source, Target, Prob)
    df <- dplyr::bind_rows(df, df0)
  }

  pcc_obj@spatial$sp.net <- list()
  pcc_obj@spatial$sp.net$prob <- adj_prob
  pcc_obj@spatial$sp.net$df <- df

  message("The communication strength based on spatial has been calculated and the result saved in the `sp.net` slot.")
  return(pcc_obj)

}

