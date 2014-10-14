#' TSCANui
#' 
#' Launch the TSCAN user interface in local machine
#'
#' This function will automatically launch the TSCAN user interface in a web browser. 
#' The user interface provides many powerful functions which is not available by command line programming.
#' It also provides a much easier and more convenient way to quickly explore single cell data and construct pseudotime analysis.
#' The user interface can also be accessed by http://zhiji.shinyapps.io/TSCAN. Neither R nor any packages are required in this online version.
#' However, it is highly recommended that the user interface be launched locally for faster running speed.
#' 
#' @export
#' @import shiny fastICA grid igraph ggplot2 plyr combinat mgcv gplots
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' \dontrun{
#'    TSCANui()
#' }


TSCANui <- function() {
      shiny::runApp(system.file("shiny",package="TSCAN"))
}
