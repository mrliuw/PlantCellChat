#' @title Plot Ligand-Receptor Pair Statistics
#'
#' @description This function generates a bar plot of ligand-receptor pair statistics for different signaling molecules.
#'
#' @param pcc_obj A PlantCellChat object.
#' @param thresh.p Numeric, specifying the p-value threshold for significance. Default is 0.05.
#' @param label.cex Numeric, specifying the text size of axis labels in the plot.
#' @param coord_flip Logical, indicating whether to flip the x and y axes for better readability of long labels. Default is `TRUE`.
#'
#' @return A bar plot displaying the count of significant interactions for each ligand-receptor pair and signaling molecule type.
#' @export
#'
PlottingLRpairStats <- function(pcc_obj, thresh.p = 0.05, label.cex = 11, coord_flip = T) {
  df <- pcc_obj@net$df

  df$CCI <- paste(df$Source, df$Target, sep = " - ")

  df_count <- df %>%
    dplyr::filter(Pvalue <= thresh.p) %>%
    dplyr::mutate(Signal = ifelse(is.na(Signal), "Unknown", Signal)) %>%
    dplyr::group_by(CCI, Signal) %>%
    dplyr::summarise(Count = n(), .groups = "drop")

  signal <- unique(df_count$Signal)

  color_map <- c("ABA" = "#4ba79b", "BR" = "#d5f1eb", "Ca2+" = "#cae9fd", "ET" = "#f4a0a0",
                 "GA" = "#dc8bb2", "IAA" = "#ab82c4", "K+" = "#5980c9", "Mg2+" = "#007089",
                 "NO3-" = "#2f4858", "SA" = "#3d5378", "Unknown" = "#2379b5", "chitin" = "#9e5394",
                 "elf18" = "#d14d81", "flg22" = "#f2545c", "JA" = "#675790")

  colors <- color_map[names(color_map) %in% signal]

  plot <- df_count %>%
    ggplot(aes(x = CCI, y = Count, fill = Signal)) +
    geom_bar(stat = "identity") +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1, color = "black", size = label.cex),
          axis.text.y = element_text(color = "black", size = label.cex),
          axis.ticks.y = element_blank()) +
    theme(legend.title = element_text(color = "black",size = 11)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(df_count$Count) / 100) * 100 + 100)) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "", fill = "Signaling molecule")

  if (coord_flip) {
    plot <- plot + coord_flip()
  }

  print(plot)
}
