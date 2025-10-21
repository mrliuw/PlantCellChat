#' @title Plot Signal Contribution bar plot
#'
#' @description
#' This function visualizes the contribution of ligand-receptor (LR) pairs or signal contributions.
#'
#' @param pcc_obj A PlantCellChat object.
#' @param ligand.type A character string specifying the type of ligand, either "lrpairs" (ligand-receptor pairs) or "signal".
#' @param key.signal A character string specifying the key signaling molecule to filter the network.
#' @param key.source A character string specifying the cell type or tissue where the ligand located.
#' @param key.target A character string specifying the cell type or tissue where the receptor located.
#' @param label.cex Numeric, scaling factor for vertex label font size.
#' @param title.cex Numeric, scaling factor for title font size.
#' @param expand.limits Numeric, Expanding the plot limits on the y-axis. Default is -10.
#' @param input.color Character. Hex color code for the bars in the plot. Default is "#4682b4".
#' @param show.title Logical, indicating whether to display the plot title.
#'
#' @return A ggplot object representing the bar plot.
#' @export

PlottingSignalContribution <- function(pcc_obj,
                                       ligand.type = "lrpairs",
                                       key.signal = NULL,
                                       key.source = NULL,
                                       key.target = NULL,
                                       top.n = 10,
                                       label.cex = 10,
                                       title.cex = 11,
                                       expand.limits = -10,
                                       input.color = "#999999",
                                       show.title = T,
                                       coord.flip = F) {
  if (!(ligand.type %in% c("lrpairs", "signal"))) {
    stop("The `ligand.type` parameter must be either 'lrpairs' or 'signal'. Please provide a valid input.")
  }

  if (ligand.type == "lrpairs") {
    df <- pcc_obj@net$df %>%
      dplyr::mutate(Interaction_name = gsub("_", " - ", Interaction_name))

    if (!is.null(key.signal)) {
      df <- df %>%
        dplyr::filter(Signal == key.signal)

      if (!(key.signal %in% df$Signal)) {
        stop("The specified `key.signal` parameter does not exist in the data.")
      }

      if (!is.null(key.source) && !is.null(key.target)) {
        if (!(key.source %in% df$Source)) {
          stop("The specified `key.source` parameter does not exist in the data.")
        }
        if (!(key.target %in% df$Target)) {
          stop("The specified `key.target` parameter does not exist in the data.")
        }
        df <- df %>%
          dplyr::filter(Source == key.source, Target == key.target)
      } else if (!is.null(key.source) && is.null(key.target)) {
        if (!(key.source %in% df$Source)) {
          stop("The specified `key.source` parameter does not exist in the data.")
        }
        df <- df %>%
          dplyr::filter(Source == key.source)
      } else if (is.null(key.source) && !is.null(key.target)) {
        stop("The `key.target` parameter should not be specified without `key.source`.")
      }

      df_filter <- df %>%
        group_by(Interaction_name) %>%
        summarise(Prob = sum(Prob)) %>%
        dplyr::mutate(Contribution = Prob / sum(Prob),
                      Interaction_name = forcats::fct_reorder(Interaction_name, Contribution, .desc = ifelse(coord.flip == T, T, F))) %>%
        dplyr::arrange(desc(Contribution)) %>%
        dplyr::filter(row_number() <= top.n)

      if (coord.flip == F) {
        plot <- df_filter %>%
          ggplot(aes(x = Contribution, y = Interaction_name)) +
          geom_bar(stat = "identity", fill = input.color, width = 0.7) +
          labs(
            x = paste(key.source, if (!is.null(key.target)) paste("-", key.target) else "", "\nRelative Contribution"),
            y = "",
            title = if (show.title) paste("Contribution of LR-pairs to", key.signal) else NULL
          ) +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = label.cex, color = "black"),
            axis.text.x = element_text(size = label.cex *1.1, color = "black"),
            axis.title.x = element_text(size = label.cex, color = "black"),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = title.cex)
          ) +
          expand_limits(x = 0, y = expand.limits)
      }

      if (coord.flip == T) {

        plot <- df_filter %>%
          ggplot(aes(x = Contribution, y = Interaction_name)) +
          geom_bar(stat = "identity", fill = input.color) +
          labs(
            x = paste(key.source, if (!is.null(key.target)) paste("-", key.target) else "", "\nRelative Contribution"),
            y = "",
            title = if (show.title) paste("Contribution of LR-pairs to", key.signal) else NULL
          ) +
          theme_bw() +
          theme(
            axis.text.y =element_text(size = label.cex * 1.1, color = "black"),
            axis.text.x = element_text(angle = 75, hjust = 1, size = label.cex, color = "black"),
            axis.title.y = element_text(size = label.cex, color = "black"),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = title.cex)
          ) +
          scale_x_continuous(expand = c(0, 0), limits = c(0,max(df_filter$Contribution) + max(df_filter$Contribution)/10)) +
          coord_flip()
      }
      print(plot)
    } else {
      stop("Please provide `key.signal` parameter to plot the contribution of LR-pairs")
    }
  }

  if (ligand.type == "signal") {
    if (!is.null(key.source)) {
      stop("The `key.source` parameter should not be specified when `ligand.type` = 'signal'.")
    }

    df <- pcc_obj@netSignal$df %>%
      dplyr::mutate(Interaction_name = sub("^[^_]*_", "", Interaction_name))

    if (!is.null(key.signal)) {
      df <- df %>%
        dplyr::filter(Signal == key.signal)

      if (!(key.signal %in% df$Signal)) {
        stop("The specified `key.signal` does not exist in the data.")
      }

      if (!is.null(key.target)) {
        if (!(key.target %in% df$Target)) {
          stop("The specified `key.target` does not exist in the data.")
        }
        df <- df %>%
          dplyr::filter(Target == key.target)
      }

      df_filter <- df %>%
        group_by(Interaction_name) %>%
        summarise(Prob = sum(Prob)) %>%
        dplyr::mutate(Contribution = Prob / sum(Prob),
                      Interaction_name = forcats::fct_reorder(Interaction_name, Contribution, .desc = ifelse(coord.flip == T, T, F))) %>%
        dplyr::arrange(desc(Contribution)) %>%
        dplyr::filter(row_number() <= top.n)

      if (coord.flip == F) {
        plot <- df_filter %>%
          ggplot(aes(x = Contribution, y = Interaction_name)) +
          geom_bar(stat = "identity", fill = input.color) +
          labs(x = paste(key.target,"Relative Contribution"), y = "",
               title = if (show.title) paste("Contribution of receptors to", key.signal) else NULL) +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = label.cex, color = "black"),
            axis.text.x = element_text(size = label.cex * 1.1, color = "black"),
            axis.title.x = element_text(size = label.cex, color = "black"),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = title.cex)
          ) +
          expand_limits(x = 0, y = expand.limits)
      }

      if (coord.flip == T) {

        plot <- df_filter %>%
          ggplot(aes(x = Contribution, y = Interaction_name)) +
          geom_bar(stat = "identity", fill = input.color) +
          labs(x = paste(key.target,"Relative Contribution"), y = "",
               title = if (show.title) paste("Contribution of receptors to", key.signal) else NULL) +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = label.cex * 1.1, color = "black"),
            axis.text.x = element_text(angle = 80, hjust = 1, size = label.cex, color = "black"),
            axis.title.y = element_text(size = label.cex, color = "black"),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size = title.cex)
          ) +
          scale_x_continuous(expand = c(0, 0), limits = c(0,max(df_filter$Contribution) + max(df_filter$Contribution)/10)) +
          coord_flip()
      }
      print(plot)
    } else {
      stop("Please provide `key.signal` parameter to plot the contribution of receptors")
    }
  }
}
