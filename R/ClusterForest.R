#' Clustering the classification trees in a forest based on similarities
#'
#' A function to get insight into a forest of classification trees by clustering the trees in a forest using Partitioning Around Medoids (PAM, Kaufman & Rousseeuw, 2009), based on similarities (see Sies & Van Mechelen, 2020).
#'
#' This function uses the Partitioning Around Medoids (PAM) algorithm, to cluster trees in a forest,
#' based on similarities between the trees. The function takes a treesimilarities object as input (see treesimilarities() or as.treesimilarities() to create a treesimilarities object)
#' and returns a clusterforest object.
#'
#' The PAM algorithm can be run for solutions with varying number of clusters, depending on the fromclus and toclus arguments.
#' Results will be returned for each solution.
#'
#' @param similarities A treesimilarities object containing a similarity matrix and the forest object on which the similarities were calculated
#' @param fromclus The lowest number of clusters for which the PAM algorithm should be run
#' @param toclus The highest number of clusters for which the PAM algorithm should be run
#' @return The function returns an object of class clusterforest, with attributes:
#' \item{medoids}{the position of the medoid trees in the forest (i.e., which element of the list of partytrees)}
#' \item{medoidtrees}{the medoid trees}
#' \item{clusters}{The cluster to which each tree in the forest is assigned}
#' \item{avgsilwidth}{The average silhouette width for each solution (see Kaufman and Rousseeuw, 2009)}
#' \item{agreement}{For each solution, the agreement between the predicted class label for each observation based on the forest as a whole, and those based on the
#' medoids only (see Sies & Van Mechelen,2020)}
#' \item{withinsim}{Within cluster similarity for each solution (see Sies & Van Mechelen, 2020)}
#' @export
#' @importFrom cluster pam
#' @importFrom graphics axis plot
#' @import MASS
#' @references \cite{Kaufman, L., & Rousseeuw, P. J. (2009). Finding groups in data: an introduction to cluster analysis (Vol. 344). John Wiley & Sons.}
#' @references \cite{Sies, A. & Van Mechelen I. (2020). C443: An R-package to see a forest for the trees. Journal of Classification.}
#' @examples
#' require(MASS)
#' require(rpart)
#'#Function to draw a bootstrap sample from a dataset
#'DrawBoots <- function(dataset, i){
#'set.seed(2394 + i)
#'Boot <- dataset[sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE),]
#'return(Boot)
#'}
#'
#'#Function to grow a tree using rpart on a dataset
#'GrowTree <- function(x,y,BootsSample, minsplit = 40, minbucket = 20, maxdepth =3){
#'  controlrpart <- rpart.control(minsplit = minsplit, minbucket = minbucket, maxdepth = maxdepth,
#'   maxsurrogate = 0, maxcompete = 0)
#'  tree <- rpart(as.formula(paste(noquote(paste(y, "~")), noquote(paste(x, collapse="+")))),
#'  data = BootsSample, control = controlrpart)
#'  return(tree)
#'}
#'
#'#Use functions to draw 15 boostrapsamples and grow a tree on each sample
#'Boots<- lapply(1:15, function(k) DrawBoots(Pima.tr ,k))
#'Trees <- lapply(1:15, function (i) GrowTree(x=c("npreg", "glu",  "bp",  "skin",
#' "bmi", "ped", "age"), y="type", Boots[[i]] ))
#'
#'#Turn the lists of dataframes and rpart trees to a forest object
#'myforest<- forest(Pima.tr,Boots,Trees)
#'
#'#Calculate similarities between trees in forest, based on similarity measure m=1.
#'Simmatrix1 <- treesimilarities(forest=myforest, m=1)
#'
#'#cluster forest
#'ClusterForest <- clusterforest(Simmatrix1, 1, 5)



clusterforest<- function(similarities, fromclus, toclus){
  simmatrix=similarities$sim
  forest=similarities$forest
  trees=forest$partytrees
  treedata=forest$treedata
  fulldata=forest$fulldata

  X= unlist(unique(lapply(1:length(forest$partytrees), function (k) attr(forest$partytrees[[k]]$terms, "term.labels"))))
  Y= unlist(lapply(1:length(forest$partytrees), function(k) colnames(forest$treedata[[k]])[!colnames(forest$treedata[[k]]) %in% X]))


  if(class(similarities) != "treesimilarities" ) {
    cat("similaritymatrix must be a treesimilarities object")}


  medoids<- list(0)
  clusters<- list(0)
  mds<- list(0)
  sil<- list(0)
  agreement<-list(0)
  sums<- list(0)

  for (i in fromclus:toclus){
    clustering <- pam(x = 1 - simmatrix, k = i, diss = TRUE)


    medoids[[i]] = clustering $ medoids
    clusters[[i]] = clustering $ clustering

    meds<- list(0)
    for(j in 1:i){
      meds[[j]] <- trees[[medoids[[i]][j]]]
    }

    mds[[i]]<-meds
    sil[[i]] <-  clustering $ silinfo $ avg.width

    pamtrees <- lapply(1:length(trees), function (i) pamtree(fulldata,treedata[[i]], Y[i],trees[[i]]))
    g<- lapply(1:length(trees), function(k) pamtrees[[k]]$predresp)
    forestpred <- sapply(1:nrow(fulldata), function (k) levels(g[[1]])[which.max( c(sum(unlist(lapply(g, `[[`, k)) ==  levels(unlist(lapply(g, `[[`, k)))[1], na.rm=TRUE),  sum(unlist(lapply(g, `[[`, k)) ==  levels(unlist(lapply(g, `[[`, k)))[2], na.rm=TRUE) )   )] )

    gmed <- g[c(medoids[[i]])]
    w <- table(clusters[[i]])
    medpred <- sapply(1:nrow(fulldata), function (k) levels(gmed[[1]])[which.max( c(sum(unlist(lapply(gmed, `[[`, k)) ==  levels(unlist(lapply(gmed, `[[`, k)))[1], na.rm=TRUE),  sum(unlist(lapply(gmed, `[[`, k)) ==  levels(unlist(lapply(gmed, `[[`, k)))[2], na.rm=TRUE) )   )] )

    agreement[[i]] <- mean(forestpred == medpred)

    sumw<- numeric(0)
    for (j in 1:i){
      sumw[j] <- sum(simmatrix[clusters[[i]]==j, medoids[[i]][j]])
    }
    sums[[i]]<- sum(sumw) / dim(simmatrix)[1]
  }


  value <- list(medoids = medoids, medoidtrees = mds, clusters=clusters, avgsilwidth=sil, agreement=agreement, withinsim=sums)
  attr(value, "class") <- "clusterforest"
  return(value)
}





