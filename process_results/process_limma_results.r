# order results by adjusted p-value
order_results_by_pvalue <- function(results) {
  ordered_results <- results[order(results$adj.P.Val), ]
  return(ordered_results)
}

# order results by log fold change
order_results_by_logFC <- function(results) {
  ordered_results <- results[order(results$logFC, decreasing = TRUE), ]
  return(ordered_results)
}

# order results by both p-value and log fold change
order_results_by_both <- function(results) {
  ordered_results <- results[order(results$adj.P.Val, decreasing = FALSE, results$logFC, decreasing = TRUE), ]
  return(ordered_results)
}


# get top 50 results by adjusted p-value
get_top_50_by_pvalue <- function(ordered_results) {
  top_50 <- head(ordered_results, 50)
  return(top_50)
}

# Function to generate volcano plot
create_volcano_plot <- function(results, p_value_threshold = 0.05, logFC_threshold = 1, title = "Volcano Plot") {
  pdf("results_plots.pdf")
  # Filter significant genes
  significant_genes <- subset(results, adj.P.Val < p_value_threshold & abs(logFC) > logFC_threshold)
  
  # Create volcano plot
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = ifelse(adj.P.Val < p_value_threshold & abs(logFC) > logFC_threshold, "Significant", "Not Significant")), size = 2) + 
    geom_text_repel(data = significant_genes, 
                    aes(label = rownames(significant_genes)), 
                    size = 3, 
                    nudge_x = 0.2, 
                    nudge_y = 0.2,
                    segment.size = 0.2,
                    max.overlaps = Inf) + 
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) + 
    geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed") + 
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") + 
    labs(title = title, x = "Log Fold Change", y = "-log10(p-value adjusted)") + 
    theme_bw() + 
    theme(legend.position = "none")
  
  return(volcano_plot)
}
