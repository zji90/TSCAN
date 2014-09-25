#' plotpseudotime
#' 
#' Plot the TSP constructed pseudotime time
#'
#' This function will plot the gene expression data after dimension reduction and link the data points with the constructed pseudotime path. It is written by plot_spanning_tree functioin in package monocle.
#' 
#' @param pseudotimedata The exact output from \code{\link{TSPpseudotime}}.
#' @param x The column of data after dimension reduction to be plotted on the horizontal axis.
#' @param y The column of data after dimension reduction to be plotted on the vertical axis.
#' @param show_tree Whether to show the links between cells connected in the minimum spanning tree.
#' @param show_cell_names Whether to draw the name of each cell in the plot.
#' @param cell_name_size The size of cell name labels if show_cell_names is TRUE.
#' @param markerexpr The gene expression used to define the size of nodes.
#' @return A ggplot2 object.
#' @export
#' @import ggplot2 plyr grid
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @references Cole Trapnell and Davide Cacchiarelli et al (2014): The dynamics and regulators of cell fate decisions are revealed by pseudo-temporal ordering of single cells. Nature Biotechnology
#' @seealso \code{\link{TSPpseudotime}} for examples
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' #Choose STAT2 gene expression as marker gene
#' STAT2expr <- log2(lpsdata["STAT2",]+1)
#' lpspseudotime <- TSPpseudotime(procdata, geneexpr = STAT2expr, dim = 2)
#' plotpseudotime(lpspseudotime, markerexpr = STAT2expr)

plotpseudotime <- function(pseudotimedata, x = 1, y = 2, show_tree = T, show_cell_names = T, cell_name_size = 3, markerexpr = NULL) {
      source_PCA_dim_1 <- source_PCA_dim_2 <- value <- sample_name <- NULL #To overcome No visible binding for global variable Note in R CMD check
      lib_info_with_pseudo <- pseudotimedata[[1]]
      TSPorder <- lib_info_with_pseudo$sample_name
      lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
      S_matrix <- pseudotimedata[[2]]
      pca_space_df <- data.frame(S_matrix[,c(x, y) ])
      colnames(pca_space_df) <- c("PCA_dim_1", "PCA_dim_2")
      pca_space_df$sample_name <- row.names(pca_space_df)
      pca_space_with_state_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
      
      edge_list <- data.frame(TSPorder[-length(TSPorder)],c(TSPorder[2:length(TSPorder)]),stringsAsFactors = F)
      colnames(edge_list) <- c("source", "target")
      edge_df <- merge(pca_space_with_state_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
      edge_df <- rename(edge_df, c(PCA_dim_1 = "source_PCA_dim_1", PCA_dim_2 = "source_PCA_dim_2"))
      edge_df <- merge(edge_df, pca_space_with_state_df[, c("sample_name", "PCA_dim_1", "PCA_dim_2")], by.x = "target", by.y = "sample_name", all = TRUE)
      edge_df <- rename(edge_df, c(PCA_dim_1 = "target_PCA_dim_1", PCA_dim_2 = "target_PCA_dim_2"))
      if (!is.null(markerexpr)) {
            markers_exprs <- data.frame(value=markerexpr)
            edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name", by.y = "row.names")
            g <- ggplot(data = edge_df, aes(x = source_PCA_dim_1, y = source_PCA_dim_2, size = log10(value + 0.1)))
      } else {
            g <- ggplot(data = edge_df, aes(x = source_PCA_dim_1, y = source_PCA_dim_2))
      }
      if (show_tree) {
            g <- g + geom_segment(aes_string(xend = "target_PCA_dim_1", yend = "target_PCA_dim_2", color = "State"), size = 0.3, linetype = "solid", na.rm = TRUE)
      }
      g <- g + geom_point(aes_string(color = "State"), na.rm = TRUE)
      if (show_cell_names) {
            g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
      }
      g <- g + theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            labs(title="Pseudotime ordering plot") + ylab("Component 1") + xlab("Component 2") + theme(legend.position = "top", 
                                                                                                       legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
            theme(panel.background = element_rect(fill = "white"))            
      g       
}
