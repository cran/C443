#' Check empirical equivalences between predictors
#'
#' Function to check for empirical equivalence relations between predictors in a data set by visualizing (partial) correaltions
#' @param data Original full data set
#' @param X The names of the predictors in the data set
#'
#' @return \item{cp}{Heatmap of correlations between all predictors} 
#' \item{Graph_pcor}{Partial correlation network between all predictors} 
#' \item{Graph_pcor_bon}{Partial correlation network between all predictors, with bonferroni correction} 
#' \item{corMat}{Matrix with correlations between all variables}
#' @export
#' @import MASS
#' @importFrom ggplot2 qplot 
#' @importFrom reshape2 melt
#' @importFrom qgraph qgraph cor_auto qgraph
#' @examples
#' require(MASS)
#' emp_eq(data = Pima.tr, X = c("npreg", "glu", "bp", "skin", "bmi", "ped", "age"))

emp_eq <- function(data, X){
  Var1 <- Var2 <- value <- NULL
  x <- data[, X]
  corMat <- cor_auto(x)
  meltedcorMat<- melt(corMat)
  cp <- qplot(x = Var1, y = Var2, data = meltedcorMat, fill = value, geom = "tile", main = "Correlations (Empirical Equivalence)", xlab = "", ylab = "")
  Graph_pcor <- qgraph(corMat, graph = "cor", layout = "spring")
  Graph_pcor_bon <- qgraph(corMat, graph = "cor", layout = "spring", threshold = "bonferroni",
                       sampleSize = nrow(data), alpha = 0.05)
  return(list(cp=cp, Graph_pcor, Graph_pcor_bon, corMat))
}
