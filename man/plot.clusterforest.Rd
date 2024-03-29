% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterMethods.R
\name{plot.clusterforest}
\alias{plot.clusterforest}
\title{Plot a clusterforest object}
\usage{
\method{plot}{clusterforest}(x, solution = NULL, predictive_plots = FALSE, ...)
}
\arguments{
\item{x}{A clusterforest object}

\item{solution}{The solution to plot the medoid trees from. If NULL, plots with the average silhouette width, within cluster similiarty
(and predictive accuracy) per solution are returned. Default = NULL}

\item{predictive_plots}{Indicating whether predictive plots should be returned: A plot showing the predictive accuracy
when making predictions based on the medoid trees, and a plot of the agreement between the class label
for each object predicted on the basis of the random forest as a whole versus based on the medoid trees. Default = FALSE.}

\item{...}{Additional arguments that can be used in generic plot function, or in plot.party.}
}
\description{
A function that can be used to plot a clusterforest object, either by returning plots
with information such as average silhouette width and within cluster siiliarity on the cluster solutions, 
or plots of the medoid trees of each solution.
}
\details{
This function can be used to plot a clusterforest object in two ways. If it's used without specifying a solution,
then the average silhouette width, and within cluster similarity measures are plotted for each solution. 
If additionally, predictive_plots=TRUE, two more plots are returned, namely a plot showing for each solution the 
predictive accuracy when making predictions based on the medoid trees, and a plot showing for each solution the agreement between
the class label for each object predicted on the basis of the random forest as a whole versus based on the medoid trees.
These plots may be helpful in deciding how many clusters are needed to summarize the forest (see Sies & Van Mechelen, 2020).

If the function is used with the clusterforest object and the number of the solution, then the medoid tree(s)
of that solution are plotted.
}
\examples{
require(MASS)
require(rpart)
#Function to draw a bootstrap sample from a dataset
DrawBoots <- function(dataset, i){
set.seed(2394 + i)
Boot <- dataset[sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE),]
return(Boot)
}

#Function to grow a tree using rpart on a dataset
GrowTree <- function(x,y,BootsSample, minsplit = 40, minbucket = 20, maxdepth =3){
 controlrpart <- rpart.control(minsplit = minsplit, minbucket = minbucket,
 maxdepth = maxdepth, maxsurrogate = 0, maxcompete = 0)
 tree <- rpart(as.formula(paste(noquote(paste(y, "~")),
 noquote(paste(x, collapse="+")))), data = BootsSample,
 control = controlrpart)
 return(tree)
}

#Use functions to draw 20 boostrapsamples and grow a tree on each sample
Boots<- lapply(1:10, function(k) DrawBoots(Pima.tr ,k))
Trees <- lapply(1:10, function (i) GrowTree(x=c("npreg", "glu",  "bp",
 "skin",  "bmi", "ped", "age"), y="type",
Boots[[i]] ))

ClusterForest<- clusterforest(observeddata=Pima.tr,treedata=Boots,trees=Trees,m=1,
fromclus=1, toclus=5, sameobs=FALSE, no_cores=2)
plot(ClusterForest)
plot(ClusterForest,2)
}
\references{
\cite{Sies, A. & Van Mechelen I. (2020). C443: An R-package to see a forest for the trees. Journal of Classification.}
}
