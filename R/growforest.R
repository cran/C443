#' Grow a forest of classification trees or tree-based treatment regimes
#'
#' Function to grow a forest based on (a) one data set and one outcome variable, by drawing bootstrap samples and growing a tree on each bootstrap sample, (b) one data set and multiple outcome variables, by drawing bootstrap samples and growing a tree for each outcome variable on each bootstrap sample, (c) multiple data sets and one outcome variable, by growing a tree on each data set. Trees can be either classification trees estimated using rpart or tree-based treatment regimes estimated using the method of Zhang et al (2012) with the AIPWE and rpart
#' @param data The data set from which bootstrap samples should be drawn, or the data sets on which the trees should be grown
#' @param X The names of the predictor variables in the data set that will be used as possible split variables
#' @param Y The name of the outcome variable(s) in the data set
#' @param ntrees The number of trees that should be grown on each data set (i.e., the number of bootstrap samples that should be drawn)
#' @param regmod NULL by default, in case of a treatment regime, it should contain the outcome model that should be used for the augmentation
#' @param A NULL by default, in case of a treatment regime, it should denote the name of the variable that indicates the assigned treatment alternative in the data set 
#' @param regime FALSE by default, TRUE if tree-based treatment regimes instead of classification trees are desired.
#' @param minsplit 40 by default, indicates the he minimum number of observations that must exist in a node in order for a split to be attempted.
#' @param minbucket 20 by default, indicates the minimum number of observations in any terminal node.
#' @param maxdepth 3 by default, the maximum depth of any node of the final tree, with the root node counted as depth 0.
#'
#' @return \item{partytrees}{The classification trees or tree-based treatment regimes saved as party objects}
#' \item{Boots}{The drawn bootstrap samples on which the trees/treatment regimes were based}
#' @export
#' @importFrom rpart rpart.control rpart prune.rpart
#' @importFrom partykit as.party
#' @importFrom stats lm predict as.formula 
#' @import MASS
#' @references \cite{Zhang, B., Tsiatis, A. A., Davidian, M., Zhang, M., & Laber, E. (2012). Estimating optimal treatment regimes from a classification perspective. Stat, 1(1), 103-114.}
#' @examples 
#' require(MASS)
#' 
#' # Create forest by drawing bootstrap samples and growing a tree on each bootstrap sample
#' forest <- growforest(data = Pima.tr, X = c("npreg", "glu",  "bp",  "skin",  "bmi", "ped", "age"), 
#' Y = "type", ntrees = 50)
#' 
#' # Create forest by drawing bootstrap samples and growing a tree for each outcome variable
#' # on each bootstrap sample
#' forest <- growforest(data= drugs, X =c ("Age", "Gender", "Edu", "Neuro", "Extr", "Open", "Agree", 
#' "Consc", "Impul","Sensat"), Y = c("Amphet", "Benzos", "Coke", "Ecst", "Leghighs", "LSD", "Mush",
#'  "Amyl", "Ket"), ntrees = 8)


growforest <- function (data, X, Y, ntrees, regmod = NULL, A = NULL, regime = FALSE, minsplit = 40, minbucket = 20, maxdepth =3 ){
  trees <- list(0)
  g <- list(0)
  partytrees <- list(0)
  Boots <- list(0)
  partytrees_ <- list(0)
  Boots_ <- list(0)
  
  
  if( ! is.data.frame(data)){
    ds <- length(data)
  }else{
    ds <- 1
   }

  ov <- length(Y)
  do <- max(ds, length(Y)) 
  
  for(k in 1:do){
    if(ds > 1){dataset <- data[[k]]
    }else{
      dataset<- data
     }
    if(ov > 1){
      YY <- Y[k]
      rg <- regmod[[k]]
    }else{
      YY<- Y
      rg<- regmod
     }
    
    for (i in 1:ntrees){
      set.seed(2394 + i)
      Boot <- dataset[sample(1:nrow(dataset), size = nrow(dataset), replace = T),]  # draw sample with replacement
      x <- Boot[, X]
      y <- Boot[, YY]
      
      
      if(regime == TRUE){
        a <- Boot[, A]
        regr<- lm(rg, data = Boot)  # estimate outcome model
        
        #Propensity model
        PR <- sum(a == 1) / length(a)  # estimate empirical proportions
        
        #Calculate mu1 and mu0
        Boot1 <- Boot
        Boot1[, A] <- 1
        Boot0 <- Boot
        Boot0[, A] <- 0
        
        mu0 <- predict(regr, as.data.frame(Boot0))  #predict outcome under treatment 0
        mu1 <- predict(regr, as.data.frame(Boot1))  # predict outcome under treatment 1
        
        # Estimate the contrast function
        Caipwe <- a * y / PR - (a - PR) / PR * mu1 - (1 - a) * y / (1 - PR) - (a - PR) / (1 - PR) * mu0   #calculate caipwe
        Z <- (sign(Caipwe > 0))
        W <- abs(Caipwe)   #Use absolute value as weights
        Z <- as.factor(Z)  #use sign as class label
        datacaipwe <- data.frame(x, Z, W)
      }
      
      if(regime == FALSE){
        Z <-  y
        W <- NULL
        datacaipwe <- data.frame(x, Z)
      }
      
      
      controlrpart <- rpart.control(minsplit = minsplit, minbucket = minbucket, maxdepth = maxdepth, maxsurrogate = 0, maxcompete = 0)
      tree <- rpart(as.formula(paste("Z~", noquote(paste(X, collapse="+")))), weights = W, data = datacaipwe, control = controlrpart) 
      
      
      if(nrow(tree $ frame) > 1){
        cptable <- tree $ cptable
        cp <- cptable[order(cptable[, 4]), ]
        maxcp <- cp[1, 4] + 2 * cp[1, 5]
        b <- cp[, 4]
        cpn <- as.numeric(which.max(1 / (maxcp - b)))
        cp <- cp[cpn, 1]
        tree1 <- prune.rpart(tree, cp = cp) 
      }else {
        tree1 <- tree
       }
      
      Boots[[i]] <- Boot
      partytree <- as.party(tree1)
      partytrees[[i]] <- partytree
      
    }
    Boots_[[k]] <- Boots
    partytrees_[[k]] <- partytrees
  }
  
  return(list(partytrees = unlist(partytrees_, recursive = FALSE), Boots = unlist(Boots_, recursive = FALSE)))
}

