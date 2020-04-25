#' Constructor of the forest class
#'
#' A function that coerces a list of party trees, trees that can be coerced to party trees, or classes inheriting from party (see Hothorn & Zeileis, 2015),
#' a list of dataframes on which these trees were grown, and the original data set underlying the forest, into an object of class forest.
#'
#' Objects of class forest contain three attributes: partytrees, treedata and fulldata. The first is a list of trees of class party (or a class inheriting from party),
#' the second is a list of dataframes on which the trees are based, and the last is a dataframe containing the original data set.
#' There are a couple of methods that can be used on a forest object, such as plot.forest(), print.forest(), summary.forest() and predict.forest().
#' @param fulldata The full/original dataset
#' @param treedata A list of dataframes on which the trees are based
#' @param trees A list of trees of class party, classes inheriting from party (e.g., glmtree), or classes that can be coerced to party (i.e., rpart, Weka_tree, XMLnode).
#' @return {The function returns an object of class forest with three attributes:}\item{partytrees}{A list of classification trees strored as party objects}
#' \item{treedata}{A list of dataframes on which the trees are based}\item{fulldata}{The full or original dataset underlying the forest}
#' @export
#' @import partykit
#' @import MASS
#' @import rpart
#' @references \cite{Hothorn, T., & Zeileis, A (2015). partykit: A Modular Toolkit for Recursive Partytioning in R. Journal of Machine Learning Research,16(118),3905-3909}
#' @examples
#' require(MASS)
#' require(rpart)
#' #Function to draw a bootstrap sample from a dataset
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
#'  data = BootsSample, control = controlrpart)
#'  return(tree)
#'}
#'
#'#Use functions to draw 20 boostrapsamples and grow a tree on each sample
#'Boots<- lapply(1:20, function(k) DrawBoots(Pima.tr ,k))
#'Trees <- lapply(1:20, function (i) GrowTree(x=c("npreg", "glu",  "bp",  "skin",
#'"bmi", "ped", "age"),y="type", Boots[[i]] ))
#'
#'
#'#Turn the lists of dataframes and rpart trees to a forest object
#'myforest<- forest(Pima.tr,Boots,Trees)
#'


forest <- function (fulldata, treedata, trees){
  if(typeof(trees) != "list" ) {
    cat("trees must be a list of party tree objects")
    return(NULL)
  }

  if(!'party' %in% class(trees[[1]])){
    if(class(trees[[1]])[1] == "rpart" ||class(trees[[1]])[1] == "Weka_tree"|| class(trees[[1]])[1] == "XMLnode"){
      trees<- lapply(1:length(trees), function (i) as.party(trees[[i]]))
    } else{cat("trees must be a list of party tree objects or objects that can be coerced to party trees")
      return(NULL)}

  }

  if(typeof(treedata) != "list" || class(treedata[[1]]) != "data.frame") {
    cat("data must be a list of data frames on which the trees were grown")
    return(NULL)
  }

  if(length(treedata) != length(trees)){
    cat("the number of data frames provided must be the same as the number of trees")
    return(NULL)
  }

  value <- list(partytrees = trees, treedata = treedata, fulldata=fulldata)
  class(value) <- append(class(value), "forest")
  return(value)
}

