#' Clustering the classification trees in a forest 
#'
#' Function to cluster classification trees in a forest using Partitioning Around Medoids (PAM, Kaufman & Rousseeuw, 2009).
#' @param simmatrix Similarity matrix containing the similarities between all pairs of trees in the forest
#' @param trees A list with all trees that should be clustered, each tree should be stored as party object
#' @param fulldata The original full dataset
#' @param treedata  A list with data sets on which the trees in the forest were based (i.e., one data set for each tree)
#' @param Y A vector with the name of the outcome variable on which each tree in the forest was based
#' @param fromclus The lowest number of clusters for which the clustering should be done
#' @param toclus The highest number of clusters for which the clustering should be done
#' @param A  by default, in case of a treatment regime, it should denote the name of the variable that indicates the assigned treatment alternative in the data set 

#'
#' @return \item{medoids}{the position of the medoid trees in the forest (i.e., which element of the list of trees)}
#' \item{mds}{the medoid trees}
#' \item{clusters}{The cluster number to which each tree is assigned}
#' \item{silplot}{Plot of the average silhouette width for each solution}
#' \item{withinplot}{Plot of the within cluster similarity for each solution}
#' \item{agreementplot}{Plot of the agreement between the assignments of the forest as a whole, and those based on the medoids for each solution}
#' @export
#' @importFrom cluster pam
#' @importFrom graphics axis plot
#' @import MASS
#' @references \cite{Kaufman, L., & Rousseeuw, P. J. (2009). Finding groups in data: an introduction to cluster analysis (Vol. 344). John Wiley & Sons.}
#' @references \cite{Ball, G. H., & Hall, D. J. (1965). ISODATA, a novel method of data analysis and pattern classification. Stanford research inst Menlo Park CA.}
#' @references \cite{Kaufman, L., & Rousseeuw, P. J. (2009). Finding groups in data: an introduction to cluster analysis (Vol. 344). John Wiley & Sons.}
#' @examples
#' require(MASS)
#' 
#' #Create a forest by drawing 10 bootstrap samples and growing a classification tree on each one
#' forest <- growforest(data = Pima.tr, X = c("npreg", "glu",  "bp",  "skin",  "bmi", "ped", "age"), 
#' Y = "type", ntrees = 10)
#' 
#' #Calculate similarities between all pairs of trees in the forest 
#' simmatrix <- similarities(fulldata = Pima.tr, treedata = forest[[2]], Y = rep("type", 10),
#' X = c("npreg", "glu", "bp", "skin", "bmi", "ped", "age"), trees=forest[[1]], m = 1, weight = 0)
#' 
#' #Cluster the trees in the forest.
#' cforest<- clusterforest(simmatrix = simmatrix, trees = forest[[1]], 
#' fulldata= Pima.tr, treedata=forest[[2]], Y=rep("type", 10), fromclus=1, toclus=5)
#' 
#' #Inspect medoids of five cluster solution
#' cforest$mds[[5]]

clusterforest <- function (simmatrix, trees, fulldata, treedata, Y, A=NULL, fromclus, toclus){
  medoids<- list(0)
  clusters<- list(0)
  mds<- list(0)
  sums<- list(NA)
  sil<- list(NA)
  agreement<- list(NA)
  
  for (i in fromclus:toclus){
    clustering <- pam(x = 1 - simmatrix, k = i, diss = TRUE)
    
    
    medoids[[i]] = clustering $ medoids
    clusters[[i]] = clustering $ clustering
    
    meds<- list(0)
    for(j in 1:i){
      meds[[j]] <- trees[[medoids[[i]][j]]]
    }
    
    mds[[i]]<-meds 
    
    # Within cluster mean similarity
     sumw<- numeric(0)
     for (j in 1:i){
      sumw[j] <- sum(simmatrix[clusters[[i]]==j, medoids[[i]][j]])
    }
    sums[[i]]<- sum(sumw) / dim(simmatrix)[1]
    
    
    # average silhouette width
    sil[[i]] <-  clustering $ silinfo $ avg.width
    
    
    #agreement medoids and forest
    
    tree <- as.pam(fulldata, treedata, Y, trees, A)
    g <- tree[[2]]
    forestpred <- sapply(1:nrow(fulldata), function (k) mean(as.numeric(lapply(g, `[[`, k))))
    forestpred[forestpred < 0.5] <- 0
    forestpred[forestpred > 0.5] <- 1
    gmed <- g[c(medoids[[i]])]
    w <- table(clusters[[i]])
    medpred <- sapply(1:nrow(fulldata), function (k) sum(as.numeric(lapply(gmed, `[[`, k)) * w) / length(trees))
    medpred[medpred < 0.5] <- 0
    medpred[medpred > 0.5] <- 1
    agreement[[i]] <- mean(forestpred == medpred) 
    }
  
  #### Silhouete plot
  sil[unlist(lapply(sil , is.null))] <- NA
  sil<- unlist(sil)
  silplot <- plot(sil, main = "Silhouette plot", xlab = "Number of clusters", ylab = "Average Silhouette width", xlim=c(1,toclus), xaxt="n")
  silplot <- silplot + axis(1, at = seq(from = 1, to = toclus, by = 1))
  
  # Within plot
  sums[unlist(lapply(sums , is.null))] <- NA
  M<- unlist(sums)
  withinplot <- plot(M, main="Within-cluster similarity", xlab="Number of clusters", ylab="Within cluster similarity", xlim=c(1,toclus), xaxt="n")
  withinplot<-withinplot + axis(1, at = seq(from = 1, to = toclus, by = 1))
  
  ## agreement
  agreement[unlist(lapply(agreement , is.null))] <- NA
  agreement<- unlist(agreement)
  agreementplot <- plot(agreement, main= "Agreement between predicted class labels by medoids and by forest as a whole", xlab = "Number of clusters", ylab = "Agreement", xlim = c(1,toclus), xaxt = "n", )
  agreementplot<- agreementplot + axis(1, at = seq(from = 1, to = toclus, by = 1))
  
  
  return(list(medoids = medoids, mds = mds, clusters = clusters, silplot=silplot, withinplot=withinplot, agreementplot=agreementplot ))
}





