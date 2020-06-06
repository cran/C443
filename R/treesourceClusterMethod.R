#' Mapping the tree clustering solution to a known source of variation underlying the forest
#'
#' A function that can be used to get insight into a clusterforest solution, in the case that there is a known
#' source of variation underlying the forest.
#' It visualizes the number of trees from each source that belong to each cluster.
#' @param clusterforest The clusterforest object
#' @param solution The solution
#' @return \item{multiplot}{For each value of the source, a bar plot with the number of trees that belong to each cluster}
#' \item{heatmap}{A heatmap with for each value of the source, the number of trees that belong to each cluster}
#' @export
#' @importFrom plyr mapvalues
#' @examples
#' require(rpart)
#' data_Amphet <-drugs[,c ("Amphet","Age", "Gender", "Edu", "Neuro", "Extr", "Open", "Agree",
#' "Consc", "Impul","Sensat")]
#' data_cocaine <-drugs[,c ("Coke","Age", "Gender", "Edu", "Neuro", "Extr", "Open", "Agree",
#'                          "Consc", "Impul","Sensat")]
#' data_LSD <- drugs[,c ("LSD","Age", "Gender", "Edu", "Neuro", "Extr", "Open", "Agree",
#'                          "Consc", "Impul","Sensat")]
#'
#'
#'#Function to draw a bootstrap sample from a dataset
#'DrawBoots <- function(dataset, i){
#'set.seed(2394 + i)
#'Boot <- dataset[sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE),]
#'return(Boot)
#'}
#'
#'#Function to grow a tree using rpart on a dataset
#'GrowTree <- function(x,y,BootsSample, minsplit = 40, minbucket = 20, maxdepth =3){
#'
#'  controlrpart <- rpart.control(minsplit = minsplit, minbucket = minbucket, maxdepth = maxdepth,
#'  maxsurrogate = 0, maxcompete = 0)
#'  tree <- rpart(as.formula(paste(noquote(paste(y, "~")), noquote(paste(x, collapse="+")))),
#'   data = BootsSample, control = controlrpart)
#'  return(tree)
#'}
#'
#' #Draw bootstrap samples and grow trees
#' BootsA<- lapply(1:5, function(k) DrawBoots(data_Amphet,k))
#' BootsC<- lapply(1:5, function(k) DrawBoots(data_cocaine,k))
#' Boots = c(BootsA,BootsC)
#'
#' TreesA <- lapply(1:5, function (i) GrowTree(x=c ("Age", "Gender", "Edu", "Neuro",
#' "Extr", "Open", "Agree","Consc", "Impul","Sensat"), y="Amphet", BootsA[[i]] ))
#' TreesC <- lapply(1:5, function (i) GrowTree(x=c ( "Age", "Gender", "Edu", "Neuro",
#' "Extr", "Open", "Agree", "Consc", "Impul","Sensat"), y="Coke", BootsC[[i]] ))
#' Trees=c(TreesA,TreesC)
#'
#'#Cluster the trees
#'ClusterForest<- clusterforest(fulldata=drugs,treedata=Boots,trees=Trees,m=1,
#'fromclus=2, toclus=2, treecov=rep(c("Amphet","Coke"),each=5))
#'
#' #Link cluster result to known source of variation
#' treesource(ClusterForest, 2)
treesource <- function(clusterforest, solution){
  UseMethod("treesource",clusterforest)
}

#' Mapping the tree clustering solution to a known source of variation underlying the forest
#'
#' A function that can be used to get insight into a clusterforest solution, in the case that there is a known
#' source of variation underlying the forest.
#' It visualizes the number of trees from each source that belong to each cluster.
#' @param clusterforest The clusterforest object
#' @param solution The solution
#' @export


treesource.default <- function(clusterforest, solution)
{
  print("Make sure that the clustering argument is an object from class clusterforest.")
}

#' Mapping the tree clustering solution to a known source of variation underlying the forest
#'
#' A function that can be used to get insight into a clusterforest solution, in the case that there is a known
#' source of variation underlying the forest.
#' It visualizes the number of trees from each source that belong to each cluster.
#' @param clusterforest The clusterforest object
#' @param solution The solution
#' @export
#' @importFrom plyr mapvalues
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn geom_bar ggtitle ylim
#' @importFrom stats frequency
#'
treesource.clusterforest <- function(clusterforest, solution)
{
  Clusters <- Sources <- freq <- cluster<- NULL
  source <- clusterforest$treecov
  treesource<- as.numeric(mapvalues(source, from=c(unique(source)), to=seq(1,length(unique(source)))))

  clustering <- clusterforest$clusters[[solution]]
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
