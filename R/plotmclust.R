#' plotmclust
#' 
#' Plot the model-based clustering results
#'
#' This function will plot the gene expression data after dimension reduction and show the clustering results.
#' 
#' @param mclustobj The exact output of \code{\link{exprmclust}} function.
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
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' lpsmclust <- exprmclust(procdata)
#' plotmclust(lpsmclust)

plotmclust <- function(mclustobj, x = 1, y = 2, show_tree = T, show_cell_names = T, cell_name_size = 3, markerexpr = NULL) {
      color_by = "State"
      
      lib_info_with_pseudo <- data.frame(State=mclustobj$clusterid,sample_name=names(mclustobj$clusterid))
      lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
      S_matrix <- mclustobj$pcareduceres
      ica_space_df <- data.frame(S_matrix[,c(x, y)])
      colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
      ica_space_df$sample_name <- row.names(ica_space_df)
      edge_df <- merge(ica_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
      
     
      if (!is.null(markerexpr)) {
            g <- ggplot(data = edge_df, aes(x = ICA_dim_1, y = ICA_dim_2, size = log10(markerexpr + 0.1)))
      } else {
            g <- ggplot(data = edge_df, aes(x = ICA_dim_1, y = ICA_dim_2))
      }
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
      if (show_cell_names) {
            g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
      }
      
      if (show_tree) {
            clucenter <- mclustobj$clucenter
            clulines <- NULL
            allsp <- shortest.paths(mclustobj$MSTtree)
            longestsp <- which(allsp == max(allsp), arr.ind = T)
            MSTorder <- get.shortest.paths(mclustobj$MSTtree,longestsp[1,1],longestsp[1,2])$vpath[[1]]       
            for (i in 1:(length(MSTorder)-1)) {
                  clulines <- rbind(clulines, c(clucenter[MSTorder[i],c(x,y)],clucenter[MSTorder[i+1],c(x,y)]))
            }
            clulines <- data.frame(x=clulines[,1],xend=clulines[,3],y=clulines[,2],yend=clulines[,4])
            g <- g + geom_segment(aes_string(x="x",xend="xend",y="y",yend="yend",size=NULL),data=clulines,size=1)
            
            clucenter <- data.frame(x=clucenter[,1],y=clucenter[,2],id=1:nrow(clucenter))
            g <- g + geom_text(aes_string(label="id",x="x",y="y",size=NULL),data=clucenter,size=10)
            
      }            
      g <- g + theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(legend.position = "top", legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
            theme(panel.background = element_rect(fill = "white"))
      g       
}
