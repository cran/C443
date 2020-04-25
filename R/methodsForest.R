#' Predict method for forest objects
#'
#' A function to make predictions for observations in a dataset based on a forest object.
#'
#' This function returns predicted class labels for each observation, obtained by predicting the classification label for an observation based on each tree in the forest,
#' and subsequently taking a majority vote.
#' @param object An object of class forest
#' @param newdata A data frame containing data for which to make predictions.
#' @param ... Additional arguments
#' @return produces a vector of predicted class labels
#' @export
#' @import partykit
#' @examples
#' require(MASS)
#' require(rpart)
#' #Function to draw a bootstrap sample from a dataset.
#'DrawBoots <- function(dataset, i){
#'set.seed(2394 + i)
#'Boot <- dataset[sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE),]
#'return(Boot)
#'}
#'
#'#Function to grow a tree using rpart on a dataset.
#'GrowTree <- function(x,y,BootsSample, minsplit = 40, minbucket = 20, maxdepth =3){
#'  controlrpart <- rpart.control(minsplit = minsplit, minbucket = minbucket,
#'  maxdepth = maxdepth, maxsurrogate = 0, maxcompete = 0)
#'  tree <- rpart(as.formula(paste(noquote(paste(y, "~")), noquote(paste(x, collapse="+")))),
#'  data = BootsSample, control = controlrpart)
#'  return(tree)
#'}
#'
#'#Use functions to draw 20 boostrapsamples and grow a tree on each sample.
#'Boots<- lapply(1:20, function(k) DrawBoots(Pima.tr ,k))
#'Trees <- lapply(1:20, function (i) GrowTree(x=c("npreg", "glu",  "bp",  "skin",
#'"bmi", "ped", "age"),y="type", Boots[[i]] ))
#'
#'
#'#Turn the lists of dataframes and rpart trees to a forest object
#'myforest<- forest(Pima.tr,Boots,Trees)
#'myforest$partytrees
#'
#'#Predict labels for observations in test set Pima.te, based on the forest
#'predict(myforest, Pima.te)
predict.forest<- function(object, newdata, ...){
  predresp<- lapply(1:length(object$partytrees), function(k) predict(object$partytrees[[k]], newdata, type="resp"))

    changepr=function(predresponse){
    predresponse[predresponse<0.5] <- 0
    predresponse[predresponse!=0] <- 1
    predresponse<-factor(predresponse, levels=c(0,1))
  }


  ## glmtree prediction gives probabilities so check whether probabilities vs classes are returned
  if(class(predresp[[1]])!= "factor"){
    predresp<- lapply(1:length(object$partytrees), function (k) changepr(predresp[[k]]) )
  }


  forestpred <- sapply(1:nrow(newdata), function (k) levels(predresp[[1]])[which.max( c(sum(unlist(lapply(predresp, `[[`, k)) ==  levels(unlist(lapply(predresp, `[[`, k)))[1], na.rm=TRUE),  sum(unlist(lapply(predresp, `[[`, k)) ==  levels(unlist(lapply(predresp, `[[`, k)))[2], na.rm=TRUE) )   )] )

  predictions <- factor(forestpred, levels=c(levels(predresp[[1]])[1], levels(predresp[[1]])[2]))
  return(predictions)
}


#' Plot a forest object
#'
#' A function to plot a tree from a forest.
#' @param x A clusterforest object
#' @param treetoplot The tree from the forest to be plotted. Default=1
#' @param ... Additional arguments
#' @export
plot.forest <- function(x, ..., treetoplot=1) {
  plot(x$partytree[[treetoplot]])
}


#' Print a forest object
#'
#' A function to print a forest, or a tree from a forest.
#' @param x A clusterforest object
#' @param treetoprint The tree to print. If NULL then all trees are printed. Default=NULL.
#' @param ... Additional arguments
#' @export
#' @importFrom stats predict
print.forest <- function(x, ..., treetoprint = NULL) {
   if (is.null(treetoprint)){
     print(x$partytree)
   } else{
   print(x$partytree[[treetoprint]])}
}


#' Summarize a forest object
#'
#' A function to summarize a forest.
#' @param object A forest object
#' @param ... Additional arguments
#' @export
summary.forest <- function(object, ...) {
 ntrees=length(object$partytrees)
 nsplits= sapply(1:ntrees, function(k) (length(object$partytrees[[k]])-1)/2)
 cat("Number of trees in the forest: ", ntrees, "\n", "Average number of splits per tree: " , mean(nsplits))
}

