#' The number of trees of each source that belong to each cluster
#'
#' Function to visualize the number of trees from each source that belong to each cluster.
#' @param source A vector with the name of the source on which each tree in the forest was based
#' @param clustering A vector with the clusternumber to which each tree belongs
#'
#' @return \item{multiplot}{For each outcome variable, a bar plot with the number of trees that belong to each cluster}
#' \item{heatmap}{A heatmap with for each outcome variable, the number of trees that belong to each cluster} 
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom gridExtra grid.arrange 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn geom_bar ggtitle ylim
#' @importFrom stats frequency
#' @examples
#' #Grow forest based on multiple outcome variables, with 5 trees for each outcome variable
#' forest <- growforest(drugs, X = c("Age", "Gender", "Edu", "Neuro", "Extr", "Open", "Agree",
#'  "Consc", "Impul","Sensat"), Y = c("Amphet", "Benzos", "Coke", "Ecst"), ntrees = 5)
#' 
#' #Calculate similarities between the trees in the forest
#' simmatrix1 <- similarities(fulldata = drugs, treedata = forest[[2]], Y = rep(c("Amphet",
#'  "Benzos", "Coke", "Ecst"), each = 5), 
#'  X = c("Age", "Gender", "Edu", "Neuro", "Extr", "Open", "Agree", "Consc", "Impul","Sensat"),
#'  trees = forest[[1]], m = 1, weight = 0)
#'
#' #Cluster the trees in the forest 
#' clusters <- clusterforest(simmatrix=simmatrix1, trees= forest[[1]], fulldata=drugs, 
#' treedata=forest[[2]], Y = rep(c("Amphet",
#' "Benzos", "Coke", "Ecst"), each = 5), 
#' fromclus=3, toclus=3)
#' 
#' #Visualize the number of trees for each source that belong to each cluster
#' treesource(source = rep(c("Amphet", "Benzos", "Coke", "Ecst"), each = 5),
#' clustering = clusters $ clusters[[3]])
 

treesource <- function(source, clustering){
  Clusters <- Sources <- freq <- cluster<- NULL
  treesource <- rep(1:length(unique(source)), each = sum(source == source[1]))
  Real <- matrix(c(0), length(unique(source)), length(unique(clustering)))
  for(i in 1:length(unique(clustering))){
    for (j in 1:length(unique(treesource))){
      Real[j,i] <- length(clustering[clustering == i & treesource == j])
    }
  }  
  row.names(Real) <- c(paste(unique(source)))
  colnames(Real) <- c(paste(unique(clustering)))
  
  # heatmap
  df.data <- expand.grid(Sources = c(paste(unique(source)))
                              , Clusters = c(paste(unique(clustering)))) 
  df.data$freq <- as.vector(Real)
  heatmap <- ggplot(data = df.data, aes(x = Clusters, y = Sources)) +
       geom_tile(aes(fill = freq))
  hm.palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space = 'Lab')  
  heatmap <- heatmap + scale_fill_gradientn(colours = hm.palette(100))
  
  ## plot per dv
  p <- list()
  for(i in 1:length(unique(source))){
    R <- data.frame(cluster = seq(1:length(unique(clustering))), frequency = Real[i, ])
    p[[i]]<- ggplot(data = R, aes(x = cluster, y = frequency)) +
      geom_bar(stat ="identity") + ylim(0,sum(as.numeric(source==unique(source[1]))))
    p[[i]] <- p[[i]] + ggtitle(paste(unique(source)[i]))
  }
  
  multiplot <- do.call(grid.arrange, p)

  return(list(multiplot = multiplot, heatmap = heatmap))
}

