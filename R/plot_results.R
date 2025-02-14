# Plot the results of the qPCR data

#' @title Plot the results of the qPCR data
#'
#' @param results A list of dataframes containing the results of the qPCR data
#'
#' @return A list of ggplot2 objects
#'
#' @import ggplot2 dplyr
#'
#' @export
#'
plot_expr <- function(results) {

  # Get the reference gene(s)
  refs <- names(results)[-1]

  # Identify the groups
  groups <- unique(results[[1]]$Group)
  if (length(groups) > 1) {
    group_by <- "Group"
  } else {
    group_by <- "Sample"
  }

  # Plot the relative expression of each sample compared to the reference sample
  plots <- list()
  for (ref in refs) {
    # Obtain the list of target genes
    targets <- unique(results[[ref]]$Target)

    # Create a plot for each target gene
    ref_plots <- list()
    for (target in targets) {
      # Filter the data
      data <- results[[ref]] %>%
        filter(Target == target)
      data[[group_by]] <- factor(data[[group_by]], levels = unique(data[[group_by]]))

      # Create the plot
      p <- ggplot(data, aes(x = !!sym(group_by), y = Relative_Expr, fill = !!sym(group_by))) +
        geom_bar(stat = "identity", position = "dodge", width = 0.8) +
        labs(title = target, y = paste0(target, "/", ref)) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )

      # Save the plot
      ref_plots[[target]] <- p
    }
    plots[[ref]] <- ref_plots
  }

  return(plots)
}
