#' Function to coerce objects to a treesimilarities object.
#'
#' A function to coerce objects to treesimilarities objects. As an argument, it takes a list with a similaritymatrix as its first element
#' and a forest object as its second.
#' This list is then coerced to a treesimilarities object (with attributes similarity matrix and forest).
#' @param obj List with as first element the similarity matrix, and as second the forest on which it was calculated
#' @param ... Additional arguments
#' @return The function returns an object of class treesimilarities. It includes two attributes:
#' forest (The forest for which similarities were calculated) and sim (similarity matrix with similarities between trees in forest based on chosen similarity measure)
#' @export
#' @import MASS
#' @import plot.matrix
#' @examples
#' require(rpart)
#' require(MASS)
#'#Function to draw a bootstrap sample from a dataset
#'DrawBoots <- function(dataset, i){
#'set.seed(2394 + i)
#'Boot <- dataset[sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE),]
#'return(Boot)
#'}
#'
#'#Function to grow a tree using rpart on a dataset
#'GrowTree <- function(x,y,BootsSample, minsplit = 40, minbucket = 20, maxdepth =3){
#'  controlrpart <- rpart.control(minsplit = minsplit, minbucket = minbucket,
#'  maxdepth = maxdepth, maxsurrogate = 0, maxcompete = 0)
#'  tree <- rpart(as.formula(paste(noquote(paste(y, "~")),
#'  noquote(paste(x, collapse="+")))), data = BootsSample, control = controlrpart)
#'  return(tree)
#'}
#'
#'#Use functions to draw 20 boostrapsamples and grow a tree on each sample
#'Boots<- lapply(1:5, function(k) DrawBoots(Pima.tr ,k))
#'Trees <- lapply(1:5, function (i) GrowTree(x=c("npreg", "glu",  "bp",  "skin",
#' "bmi", "ped", "age"), y="type", Boots[[i]] ))
#'
#'#Turn the lists of dataframes and rpart trees to a forest object
#'myforest<- forest(Pima.tr, Boots,Trees)
#'
#'#Create similarity matrix based on some measure
#'simmatrix<-matrix(c(1,0.5,0.8,0.9,0.3, 0.5,1, 0.6,0.7,0.2, 0.8,0.6,1,0.7,0.3,
#' 0.9,0.7,0.7,1,0.3, 0.3,0.2,0.3,0.3,1  ), 5,5,)
#'
#'
#'simobject<- as.treesimilarities(list(simmatrix,myforest))
as.treesimilarities <- function(obj, ...)
  UseMethod("as.treesimilarities")

#' A function to coerce objects to a treesimilarities object.
#'
#' A function to coerce objects to treesimilarities objects. It takes a list with a similaritymatrix as its first element
#' and a forest object as its second, as arguments.
#' This list is then coerced to a treesimilarities object (with attributes similarity matrix and forest).
#' @param obj List with as first element the similarity matrix, and as second the forest on which it was calculated
#' @param ... Additional arguments
#' @export
as.treesimilarities.list <- function(obj, ...) {
  sim= obj[[1]]
  forest= obj[[2]]
  value <- list(forest = forest, sim = sim)
  attr(value, "class") <- "treesimilarities"
  return(value)

}


#' print a treesimilarities object
#'
#' A function to print a treesimilarities object.
#' @param x A treesimilarities object
#' @param ... Additional arguments
#' @export
print.treesimilarities<- function(x, ...){
  print(x$sim)
}

#' Summarize a treesimilarities object
#'
#' A function to summarize a treesimilarities object.
#' @param object A treesimilarities object
#' @param ... Additional arguments
#' @export
summary.treesimilarities<- function (object, ...){
  cat("Number of trees:", ncol(object$sim), "\n")
  cat("Average Similarity per tree with all other trees: \n" , (rowSums(object$sim)-1)/nrow(object$sim))
}



#' Plot a treesimilarities object
#'
#' A function to plot a treesimilarities object.
#' @param x A treesimilarities object
#' @param ... Additional arguments
#' @export
#' @import plot.matrix
plot.treesimilarities<- function(x, ...){
  plot(x$sim, main="Tree similarity matrix")
}





