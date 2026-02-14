#' @title Merge Expression Matrices
#'
#' @description This function merges expression matrices by filling
#'     missing values with non-missing values from other matrices.
#'
#' @param matrix1 The first expression matrix.
#' @param matrix2 The second expression matrix.
#' @param matrix_list A list of expression matrices to be merged.
#'
#' @return A merged expression matrix.
#' @export
#'
#' @examples
#' matrix1 <- matrix(c(NA,2,4,NA,7,9,NA,4,11),nrow = 3)
#' matrix2 <- matrix(c(5,NA,6,21,NA,9,8,NA,5),nrow = 3)
#' matrix3 <- matrix(c(1,3,NA,7,4,NA,6,11,NA),nrow = 3)
#'
#' merged_matrix <- MergeExpMatrix(matrix1,matrix2)
#' merged_matrix <- MergeExpMatrix(matrix_list = list(matrix1,matrix2,matrix3))
#'
#' @seealso \code{\link[stats]{complete.cases}}
MergeExpMatrix <- function(matrix1 = NULL,matrix2 = NULL,matrix_list = NULL) {
  if (is.null(matrix_list)) {
    matrix_list <- list(matrix1,matrix2)
  }
  result_matrix <- matrix(NA,nrow = nrow(matrix_list[[1]]),ncol = ncol(matrix_list[[1]]))

  for (matrix in matrix_list) {
    result_matrix[which(complete.cases(matrix)),] <- matrix[which(complete.cases(matrix)),]
  }

  return(result_matrix)
}

#' @title Extract Expression Matrices
#'
#' @description This function extracts the expression matrix for specific components from
#'     the PlantCellChat object. It allows replacement of NA values or zero values with a specified value.
#'
#' @param pcc_obj A PlantCellChat object.
#' @param comps The component (or column) of the interaction data to extract expression for.
#' @param interaction.type A string defining the type of interaction. Accepts either "lri" for
#'     ligand-receptor interactions or "sri" for signaling molecule interactions.
#' @param na.replace Logical, if TRUE, replace NA values in the matrix with 1. Default is FALSE.
#' @param zero.replace Logical, if TRUE, replace NA values in the matrix with 0. Default is FALSE.
#' @param spatial Logical, if TRUE, extract expression data from the `data.signaling` slot. Default is FALSE.
#'
#' @return A matrix of expression values for the specified components, with optional replacements.
#' @export
ExtractExpMatrix <- function(pcc_obj,comps,interaction.type,na.replace = F,zero.replace = F,spatial = F) {
  interaction <- pcc_obj@lrpairs$lr.signaling
  if (interaction.type == "lri") {
    interaction <- interaction[interaction$Interaction_type == "Interacting protein",]
  }
  if (interaction.type == "sri") {
    interaction <- interaction[interaction$Interaction_type == "Signaling molecule",]
  }

  if (spatial == T) {
    exp <- as.matrix(pcc_obj@data.signaling)
  } else {
    exp <- pcc_obj@data.average
  }

  index <- interaction[,comps]
  matchrow <- match(index,rownames(exp))
  comps_matrix <- exp[matchrow,]
  comps_matrix <- matrix(as.numeric(comps_matrix),nrow = nrow(interaction),ncol = ncol(exp))
  if (na.replace == T) {
    comps_matrix[is.na(comps_matrix)] <- 1
  }
  if (zero.replace == T) {
    comps_matrix[is.na(comps_matrix)] <- 0
  }
  return(comps_matrix)
}

#' @title Correct Matrix
#'
#' @description This function corrects a matrix by replacing NA values
#'     with 1 and setting non-1 values to NA.
#'
#' @param matrix A matrix to be corrected.
#'
#' @return A corrected matrix.
#' @export
#'
#' @examples
#' set.seed(123)
#' matrix <- matrix(runif(11 * 5), nrow = 11, ncol = 5)
#' matrix[5:7,] <- NA
#' corrected_matrix <- correct(matrix)

correct <- function(matrix) {
  if (!is.matrix(matrix)) {
    stop("Input is not a matrix!")
  }
  matrix[is.na(matrix)] <- 1
  matrix[matrix != 1] <- NA
  return(matrix)
}

#' @title sub-function of `CalculateCommunStrength`
#'
#' @param pcc_obj A `PlantCellChat` object which contains the average expression data (stored in the `data.average` slot)
#'     and over-expressed ligand-receptor interaction data (stored in the `lrpairs` slot).
#' @param Kh Numeric, the dissociation constant in the Hill equation, controlling the sensitivity of ligand-receptor binding.
#' @param n Numeric, the Hill coefficient describing the cooperativity of ligand-receptor binding.
#'
#' @return The `PlantCellChat` object with an updated `prob` matrix, stored in the `net` slot, which contains the calculated
#'     communication strength probabilities for each ligand-receptor pair across cell clusters.
#'
CalculateCommunStrength.Sub <- function(pcc_obj, Kh = 0.5, n = 1) {
  if (is.null(pcc_obj@lrpairs$lr.signaling)) {
    stop("Over-expressed interaction were not extracted, please run the function 'ExtractOverExpressedInteractions' first.")
  } else {
    interaction <- pcc_obj@lrpairs$lr.signaling
    interaction <- interaction[interaction$Interaction_type == "Interacting protein",]
    if (is.null(pcc_obj@data.average)) {
      stop("Average Expression were not calculated, please run the function 'CalculateAvgExp' first.")
    } else {
      exp <- pcc_obj@data.average
    }
  }

  # Ligand
  exp_L <- ExtractExpMatrix(pcc_obj,comps = "Ligand",interaction.type = "lri",na.replace = T)

  # Receptor
  exp_single_R <- ExtractExpMatrix(pcc_obj,comps = "Receptor",interaction.type = "lri",na.replace = T)

  exp_complex1 <- ExtractExpMatrix(pcc_obj,comps = "Complex1",interaction.type = "lri")
  exp_complex2 <- ExtractExpMatrix(pcc_obj,comps = "Complex2",interaction.type = "lri")
  sub1unit <- exp_single_R * correct(exp_complex1)
  sub2unit <- (exp_single_R * exp_complex1 * correct(exp_complex2))^(1/2)
  sub3unit <- (exp_single_R * exp_complex1 * exp_complex2)^(1/3)

  exp_coI <- ExtractExpMatrix(pcc_obj,comps = "CoAreceptor",interaction.type = "lri",zero.replace = T)
  exp_coA <- ExtractExpMatrix(pcc_obj,comps = "CoIreceptor",interaction.type = "lri",zero.replace = T)

  exp_R <- MergeExpMatrix(matrix_list = list(sub1unit,sub2unit,sub3unit)) * (1 + exp_coA) * (1/(1 + exp_coI))

  # Agonist
  exp_ag <- ExtractExpMatrix(pcc_obj,comps = "Agonist",interaction.type = "lri",zero.replace = T)
  p_ag <- 1 + (exp_ag^n/(Kh^n + exp_ag^n))

  # Antagonist
  exp_anta <- ExtractExpMatrix(pcc_obj,comps = "Antagonist",interaction.type = "lri",zero.replace = T)
  p_anta <- Kh^n / (Kh^n + exp_anta^n)

  # calculate communication strength
  nLR <- nrow(interaction)
  numCluster <- nlevels(pcc_obj@idents)
  prob <- array(NA,dim = c(numCluster,numCluster,nLR))

  for (i in 1:nLR){
    prob[,,i] <- (exp_L[i,] %*% t(exp_R[i,]))^n/(Kh^n + (exp_L[i,] %*% t(exp_R[i,]))^n) * (p_ag[i,] %*% t(p_ag[i,])) * (p_anta[i,] %*% t(p_anta[i,]))
  }

  dimnames(prob)[[1]] <- levels(PccIdents(pcc_obj))
  dimnames(prob)[[2]] <- levels(PccIdents(pcc_obj))
  dimnames(prob)[[3]] <- rownames(interaction)

  pcc_obj@net$prob <- prob

  return(pcc_obj)
}

#' @title Calculate Communication Strength
#'
#' @description This function using the Law of mass to calculate the communication strength based on
#'     ligands and receptors.
#'
#' @param pcc_obj A `PlantCellChat` object which contains the average expression data (stored in the `data.average` slot)
#'     and over-expressed ligand-receptor interaction data (stored in the `lrpairs` slot).
#' @param Kh Numeric, the dissociation constant in the Hill equation, controlling the sensitivity of ligand-receptor binding.
#' @param n Numeric, the Hill coefficient describing the cooperativity of ligand-receptor binding.
#' @param exp.methods Character string, specifying the method for calculating average expression, default is "average".
#' @param num.permutations Integer, the number of permutations to run in calculating the p-value for communication strength.
#' @param seed Integer, the seed for random number generation, ensuring reproducibility.
#'
#' @return The `PlantCellChat` object with updated communication strength data stored in the `net` slot, including:
#'     - `prob`: A 3-dimensional matrix containing calculated communication probabilities for each ligand-receptor pair across clusters.
#'     - `pvalue`: A 3-dimensional matrix with p-values based on permutation tests for assessing the significance of each probability.
#'     - `df`: A data frame summarizing the calculated probabilities and p-values for each ligand-receptor interaction.
#' @export
CalculateCommunStrength <- function(pcc_obj, Kh = 0.5, n = 1, exp.methods = "average", num.permutations = 100, seed = 123) {
  if (is.null(pcc_obj@lrpairs$lr.signaling)) {
    stop("Over-expressed interaction were not extracted, please run the function 'ExtractOverExpressedInteractions' first.")
  } else {
    interaction <- pcc_obj@lrpairs$lr.signaling
    interaction <- interaction[interaction$Interaction_type == "Interacting protein",]
    if (is.null(pcc_obj@data.average)) {
      stop("Average Expression were not calculated, please run the function 'CalculateAvgExp' first.")
    } else {
      exp <- pcc_obj@data.average
    }
  }

  pcc_obj <- CalculateCommunStrength.Sub(pcc_obj, Kh, n)

  p_list <- list()

  set.seed(seed)

  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = num.permutations,
                                   complete = "=",
                                   incomplete = "-",
                                   current = ">",
                                   clear = FALSE,
                                   width = 100)

  for (i in 1:num.permutations) {
    pb$tick()
    original_label <- PccIdents(pcc_obj)
    resample_label <- sample(original_label,
                             size = length(original_label),
                             replace = F)
    suppressMessages(pcc_obj_resample <- SetPccIdent(pcc_obj,resample_label))
    suppressMessages(pcc_obj_resample <- CalculateAvgExp(pcc_obj_resample,methods = exp.methods))
    suppressMessages(pcc_obj_resample <- CalculateCommunStrength.Sub(pcc_obj_resample,Kh = Kh,n = n))

    p_list[[i]] <- pcc_obj_resample@net$prob < pcc_obj@net$prob
    p_list[[i]][p_list[[i]] == TRUE] <- 1
  }

  pcc_obj@net$pvalue <- Reduce(`+`,p_list) / num.permutations

  prob <- pcc_obj@net$prob
  pvalue <- pcc_obj@net$pvalue

  df <- data.frame()
  for (i in 1:dim(prob)[[3]]) {
    df0 <- prob[,,i] %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Target") %>%
      tidyr::gather(key = "Source",value = "Prob",-Target) %>%
      dplyr::mutate(pvalue[,,i] %>%
                      as.data.frame() %>%
                      tidyr::gather(key = "Source",value = "Pvalue") %>%
                      dplyr::select("Pvalue")) %>%
      dplyr::mutate(Interaction_name = dimnames(prob)[[3]][i]) %>%
      dplyr::left_join(interaction,by = "Interaction_name") %>%
      dplyr::select(Interaction_name, Signal, Source, Target, Prob, Pvalue)


    df <- dplyr::bind_rows(df,df0)
  }

  pcc_obj@net$df <- df

  message("The communication strength has been calculated and the result saved in the `net` slot.")

  return(pcc_obj)
}

#' @title sub-function of `CalculateSignalingStrength`
#'
#' @param pcc_obj A `PlantCellChat` object containing data on average expression (in the `data.average` slot) and
#'     over-expressed signaling interactions (in the `lrpairs` slot).
#' @param Kh Numeric, the dissociation constant in the Hill equation, which modulates the sensitivity of signaling molecule binding.
#' @param n Numeric, the Hill coefficient, which describes the cooperativity of signaling molecule binding.
#'
#' @return Returns the modified `PlantCellChat` object with an updated `prob` matrix stored in the `net` slot. This matrix
#'     contains the calculated communication strength probabilities for each ligand-receptor pair across cell clusters.
#' @export
#'
CalculateSignalingStrength.Sub <- function(pcc_obj, Kh = 0.5, n = 1, Ls = 1) {
  if (is.null(pcc_obj@lrpairs$lr.signaling)) {
    stop("Over-expressed interaction were not extracted, please run the function 'ExtractOverExpressedInteractions' first.")
  } else {
    interaction <- pcc_obj@lrpairs$lr.signaling
    interaction <- interaction[interaction$Interaction_type == "Signaling molecule",]
    if (is.null(pcc_obj@data.average)) {
      stop("Average Expression were not calculated, please run the function 'CalculateAvgExp' first.")
    } else {
      exp <- pcc_obj@data.average
    }
  }
  
  # Ligand
  # exp_L <- ExtractExpMatrix(pcc_obj,comps = "Ligand",interaction.type = "sri",na.replace = T)
  
  # Receptor
  exp_single_R <- ExtractExpMatrix(pcc_obj,comps = "Receptor",interaction.type = "sri",na.replace = T)
  
  exp_complex1 <- ExtractExpMatrix(pcc_obj,comps = "Complex1",interaction.type = "sri")
  exp_complex2 <- ExtractExpMatrix(pcc_obj,comps = "Complex2",interaction.type = "sri")
  sub1unit <- exp_single_R * correct(exp_complex1)
  sub2unit <- (exp_single_R * exp_complex1 * correct(exp_complex2))^(1/2)
  sub3unit <- (exp_single_R * exp_complex1 * exp_complex2)^(1/3)
  
  exp_coI <- ExtractExpMatrix(pcc_obj,comps = "CoAreceptor",interaction.type = "sri",zero.replace = T)
  exp_coA <- ExtractExpMatrix(pcc_obj,comps = "CoIreceptor",interaction.type = "sri",zero.replace = T)
  
  exp_R <- MergeExpMatrix(matrix_list = list(sub1unit,sub2unit,sub3unit)) * (1 + exp_coA) * (1/(1 + exp_coI))
  
  # Agonist
  exp_ag <- ExtractExpMatrix(pcc_obj,comps = "Agonist",interaction.type = "sri",zero.replace = T)
  p_ag <- 1 + (exp_ag^n/(Kh^n + exp_ag^n))
  
  # Antagonist
  exp_anta <- ExtractExpMatrix(pcc_obj,comps = "Antagonist",interaction.type = "sri",zero.replace = T)
  p_anta <- Kh^n / (Kh^n + exp_anta^n)
  
  # calculate communication strength
  nLR <- nrow(interaction)
  numCluster <- nlevels(pcc_obj@idents)
  prob <- array(NA,dim = c(1,numCluster,nLR))
  
  for (i in 1:nLR){
    prob[,,i] <- (Ls * t(exp_R[i,]))^n/(Kh^n + (Ls * t(exp_R[i,]))^n) * t(p_ag[i,]) * t(p_anta[i,])
  }
  
  dimnames(prob)[[2]] <- levels(PccIdents(pcc_obj))
  dimnames(prob)[[3]] <- rownames(interaction)
  pcc_obj@netSignal$prob <- prob
  
  return(pcc_obj)
}

#' @title Calculate Signaling Strength
#'
#' @description This function calculates the signaling strength between cell populations based on signaling molecule interactions.
#'     It uses the Law of mass to model cooperative binding effects and performs permutation testing to assess the statistical
#'         significance of each signaling strength value. The function outputs probability and p-value matrices for signaling interactions.
#'
#' @param pcc_obj A `PlantCellChat` object containing data on average expression (in the `data.average` slot) and
#'     over-expressed signaling interactions (in the `lrpairs` slot).
#' @param Kh Numeric, the dissociation constant in the Hill equation, which modulates the sensitivity of signaling molecule binding.
#' @param n Numeric, the Hill coefficient, which describes the cooperativity of signaling molecule binding.
#' @param exp.methods Character, the method for calculating average expression; default is "average".
#' @param num.permutations Integer, the number of permutations used in p-value calculation, providing statistical significance of the signaling strength.
#' @param seed Integer, a seed for random number generation, ensuring reproducibility.
#'
#' @return The `PlantCellChat` object with updated signaling strength data stored in the `netSignal` slot, including:
#'     - `prob`: A 3-dimensional matrix containing calculated signaling probabilities.
#'     - `pvalue`: A 3-dimensional matrix with p-values based on permutation tests.
#'     - `df`: A data frame summarizing the calculated probabilities and p-values for each signaling interaction.
#'
#' @export
#'
CalculateSignalingStrength <- function(pcc_obj, Kh = 0.5, n = 1, exp.methods = "average", num.permutations = 100, seed = 123, Ls = 1) {
  if (is.null(pcc_obj@lrpairs$lr.signaling)) {
    stop("Over-expressed interaction were not extracted, please run the function 'ExtractOverExpressedInteractions' first.")
  } else {
    interaction <- pcc_obj@lrpairs$lr.signaling
    interaction <- interaction[interaction$Interaction_type == "Signaling molecule",]
    if (is.null(pcc_obj@data.average)) {
      stop("Average Expression were not calculated, please run the function 'CalculateAvgExp' first.")
    } else {
      exp <- pcc_obj@data.average
    }
  }
  
  pcc_obj <- CalculateSignalingStrength.Sub(pcc_obj, Kh, n, Ls)
  
  p_list_signal <- list()
  
  set.seed(seed)
  
  pb <- progress::progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = num.permutations,
    complete = "=",
    incomplete = "-",
    current = ">",
    clear = FALSE,
    width = 100
  )
  
  for (i in 1:num.permutations) {
    pb$tick()
    
    original_label <- PccIdents(pcc_obj)
    resample_label <- sample(original_label, size = length(original_label), replace = FALSE)
    suppressMessages(pcc_obj_resample <- SetPccIdent(pcc_obj, resample_label))
    suppressMessages(pcc_obj_resample <- CalculateAvgExp(pcc_obj_resample,methods = exp.methods))
    suppressMessages(pcc_obj_resample <- CalculateSignalingStrength.Sub(pcc_obj_resample, Kh, n, Ls))
    
    p_list_signal[[i]] <- pcc_obj_resample@netSignal$prob < pcc_obj@netSignal$prob
    p_list_signal[[i]][p_list_signal[[i]] == TRUE] <- 1
  }
  
  pcc_obj@netSignal$pvalue <- Reduce(`+`, p_list_signal) / num.permutations
  
  prob <- pcc_obj@netSignal$prob
  pvalue <- pcc_obj@netSignal$pvalue
  
  df <- data.frame()
  for (i in 1:dim(prob)[[3]]) {
    source_name <- strsplit(dimnames(prob)[[3]][i], "_")[[1]][1]
    df0 <- as.data.frame(prob[,,i]) %>%
      tibble::rownames_to_column(var = "Target") %>%
      tidyr::pivot_longer(cols = -Target, names_to = "Source", values_to = "Prob") %>%
      dplyr::mutate(Source = source_name,
                    Pvalue = as.vector(pvalue[,,i]),
                    Interaction_name = dimnames(prob)[[3]][i]) %>%
      dplyr::left_join(interaction,by = "Interaction_name") %>%
      dplyr::select(Interaction_name, Signal, Source, Target, Prob, Pvalue)
    df <- dplyr::bind_rows(df, df0)
  }
  
  pcc_obj@netSignal$df <- df
  
  message("The signaling strength has been calculated and the result saved in the `netSignal` slot.")
  
  return(pcc_obj)
}
