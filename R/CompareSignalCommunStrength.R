#' @title Compare Signal Communication Strength
#'
#' @description
#' This function compares the communication strength of different groups based on the ligand-receptor pairs or signaling molecules.
#'
#'
#' @param pcc_obj A PlantCellChat object.
#' @param pcc_obj_list A list of PlantCellChat objects for comparison across multiple groups
#' @param ligand.type A character string specifying the type of ligand, either "lrpairs" (ligand-receptor pairs) or "signal".
#' @param key.signal A character string specifying the key signaling molecule to filter the network.
#' @param key.source A character string specifying the cell type or tissue where the ligand located.
#' @param key.target A character string specifying the cell type or tissue where the receptor located.
#' @param label.cex Numeric, scaling factor for vertex label font size.
#' @param title.cex Numeric, scaling factor for title font size.
#' @param expand.limits Numeric, Expanding the plot limits on the y-axis. Default is 0.
#' @param input.color Character. Hex color code for the bars in the plot.
#' @param show.title Logical, indicating whether to display the plot title.
#'
#' @return A ggplot object representing the bar plot.
#' @export
#'
CompareSignalCommunStrength <- function(pcc_obj,
                                        pcc_obj_list = NULL,
                                        ligand.type = "lrpairs",
                                        key.signal = NULL,
                                        key.source = NULL,
                                        key.target = NULL,
                                        top.n = 10,
                                        label.cex = 10,
                                        title.cex = 11,
                                        expand.limits = 0,
                                        input.color = NULL,
                                        show.title = T,
                                        coord.flip = F,
                                        scale.factor = 1) {
  if (!(ligand.type %in% c("lrpairs", "signal"))) {
    stop("The `ligand.type` parameter must be either 'lrpairs' or 'signal'. Please provide a valid input.")
  }

  default_colors <- c("#458bc9", "#ee7d79", "#FCC910", "#7CC910", "#10C9CC", "#C910CC")

  if (ligand.type == "lrpairs") {
    if (!is.null(pcc_obj_list) && is.list(pcc_obj_list)) {
      combined_df <- NULL
      for (group_name in names(pcc_obj_list)) {
        current_df <- pcc_obj_list[[group_name]]@net$df
        current_df$Group <- group_name
        combined_df <- rbind(combined_df, current_df)
      }
    }

    df <- pcc_obj@net$df %>%
      dplyr::mutate(Interaction_name = gsub("_", " - ", Interaction_name))

    combined_df <- combined_df %>%
      dplyr::mutate(Interaction_name = gsub("_", " - ", Interaction_name))

    if (!is.null(key.signal)) {
      df <- df %>%
        dplyr::filter(Signal == key.signal)

      combined_df <- combined_df %>%
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
        combined_df <- combined_df %>%
          dplyr::filter(Source == key.source, Target == key.target)

      } else if (!is.null(key.source) && is.null(key.target)) {
        if (!(key.source %in% df$Source)) {
          stop("The specified `key.source` parameter does not exist in the data.")
        }
        df <- df %>%
          dplyr::filter(Source == key.source)
        combined_df <- combined_df %>%
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

      combined_df_filter <- combined_df %>%
        filter(Interaction_name %in% df_filter$Interaction_name) %>%
        group_by(Group,Interaction_name) %>%
        summarise(Prob = sum(Prob)) %>%
        left_join(df_filter, by = "Interaction_name") %>%
        select(-Prob.y) %>%
        rename(Prob = Prob.x)

      # temp_row <- combined_df_filter[1,]
      # temp_row$Group <- "Temp"
      # temp_row$Prob <- 0
      #
      # combined_df_filter <- rbind(combined_df_filter, temp_row)

      combined_df_filter$Interaction_name <- factor(combined_df_filter$Interaction_name, levels = rev(df_filter$Interaction_name))

      combined_df_filter$Contribution <- -combined_df_filter$Contribution * scale.factor

      combined_df_filter %>%
      ggplot(aes(x = Interaction_name)) +
        geom_bar(aes(y = Prob, fill = Group), stat = "identity", width = 0.8,position = "dodge") +
        geom_bar(aes(y = Contribution), stat = "identity", fill = "#cccccc", width = 0.8) +
        coord_flip() +
        labs(
          title = if (show.title) paste(key.signal, "related LR-pairs") else NULL
        ) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = label.cex, color = "black"),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = title.cex)
        ) +
        coord_flip() +
        scale_fill_manual(
          values = setNames(c(input.color %||% default_colors)[1:length(pcc_obj_list)], names(pcc_obj_list)),
          labels = setNames(
            c(sapply(names(pcc_obj_list), function(name) paste(name, "group communication strength")),
              paste(key.source, if (!is.null(key.target)) paste("-", key.target) else "", "relative contribution")),
            names(pcc_obj_list)[1:length(pcc_obj_list)]
          )
        )
    } else {
      stop("Please provide `key.signal` parameter to plot the contribution of LR-pairs")
    }
  }

  if (ligand.type == "signal") {
    if (!is.null(pcc_obj_list) && is.list(pcc_obj_list)) {
      combined_df <- NULL
      for (group_name in names(pcc_obj_list)) {
        current_df <- pcc_obj_list[[group_name]]@netSignal$df
        current_df$Group <- group_name
        combined_df <- rbind(combined_df, current_df)
      }
    }

    if (!is.null(key.source)) {
      stop("The `key.source` parameter should not be specified when `ligand.type` = 'signal'.")
    }

    df <- pcc_obj@netSignal$df %>%
      dplyr::mutate(Interaction_name = sub("^[^_]*_", "", Interaction_name))

    combined_df <- combined_df %>%
      dplyr::mutate(Interaction_name = sub("^[^_]*_", "", Interaction_name))

    if (!is.null(key.signal)) {
      df <- df %>%
        dplyr::filter(Signal == key.signal)

      combined_df <- combined_df %>%
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
        combined_df <- combined_df %>%
          dplyr::filter(Target == key.target)
      }

      df_filter <- df %>%
        group_by(Interaction_name) %>%
        summarise(Prob = sum(Prob)) %>%
        dplyr::mutate(Contribution = Prob / sum(Prob),
                      Interaction_name = forcats::fct_reorder(Interaction_name, Contribution, .desc = ifelse(coord.flip == T, T, F))) %>%
        dplyr::arrange(desc(Contribution)) %>%
        dplyr::filter(row_number() <= top.n)

      combined_df_filter <- combined_df %>%
        filter(Interaction_name %in% df_filter$Interaction_name) %>%
        group_by(Group,Interaction_name) %>%
        summarise(Prob = sum(Prob)) %>%
        left_join(df_filter, by = "Interaction_name") %>%
        select(-Prob.y) %>%
        rename(Prob = Prob.x)

      combined_df_filter$Interaction_name <- factor(combined_df_filter$Interaction_name, levels = rev(df_filter$Interaction_name))

      combined_df_filter$Contribution <- -combined_df_filter$Contribution * scale.factor

      combined_df_filter %>%
        ggplot(aes(x = Interaction_name)) +
        geom_bar(aes(y = Prob, fill = Group), stat = "identity", width = 0.6,position = "dodge") +
        geom_bar(aes(y = Contribution), stat = "identity", fill = "#cccccc", width = 0.6) +
        coord_flip() +
        labs(
          title = if (show.title) paste(key.signal, "related Receptors") else NULL
        ) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = label.cex, color = "black"),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = title.cex)
        ) +
        coord_flip() +
        scale_fill_manual(
          values = setNames(c(input.color %||% default_colors)[1:length(pcc_obj_list)], c(names(pcc_obj_list))),
          labels = setNames(
            c(sapply(names(pcc_obj_list), function(name) paste(name, "group communication strength")),
              paste(key.target,"Relative Contribution")),
            c(names(pcc_obj_list), "Temp")[1:(length(pcc_obj_list) + 1)]
          )
        )
    } else {
      stop("Please provide `key.signal` parameter to plot the contribution of receptors")
    }
  }
}


