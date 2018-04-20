#' Calculating similarities between classification trees
#'
#' Function to calculate similarities between classification trees, based on 6 different  possible similarity measures.
#' @param fulldata The original full data set 
#' @param treedata A list with data sets on which the trees in the forest were based (i.e., one data set for each tree)
#' @param Y A vector with the name of the outcome variable on which each tree in the forest was based
#' @param X The names of the predictor variables that were used as possible split variables  
#' @param trees A list with all trees between which similarities should be computed, each tree should be stored as party object
#' @param m Similarity measure that should be used to calculate similarities, where m=1 is based on counting equal predictors or predictor-split point combinations (Equation 5 or 8  in Sies & Van Mechelen (Submitted), m=2
#' is the measure of Shannon & Banks (1999), based on counting the number of equal paths from rootnode to leafs (See Sies & Van Mechelen Submitted, Equation 2), m=3 is 
#' based on the agreement in classification lables (Chipman, 1998), see Sies & Van Mechelen (submitted), Equation 14, m=4
#' is based on the agreement of partitions (Chipman, 1998), see Sies & Van Mechelen (Submitted), Equation 13, and m=5 is based on counting equal
#' elementary conjunctions of trees transformed to disjunctive normal form (only for binary predictors, see Sies & Van Mechelen, Submitted, Equation 16). Finally M6 is based on comparing sets of predictor split piont 
#' combinations (taking into account directions of the splits) associated with a leaf, taking into account the classification label of that leaf, see Sies & Van Mechelen (submitted).
#' @param weight Indicating whether or not splitpoints should be taken into account for m=1, where 0 means no (Equation 4 in Sies & Van Mechelen, submitted) and 1 means yes (Equation 8 in Sies & Van Mechelen, submitted). 
#' @param A The name of the treatment variable in case of a forest of tree-based treatment regimes, otherwise NULL by default.
#' @param regime Indicating whether the trees in the forest are treatment regimes (TRUE) or decision trees (FALSE). Default=FALSE
#' @param tol In case that weight = 1: A vector with for each predictor the tolerance zone within which two split points of the predictor in question are assumed equal. Default=NULL
#' @return {sim}{Similarity matrix based on chosen similarity measure}
#' @export 
#' @import MASS
#' @importFrom parallel detectCores makeCluster clusterExport parSapply stopCluster
#' @importFrom igraph graph_from_incidence_matrix max_bipartite_match
#' @importFrom partykit nodeids data_party id_node kids_node varid_split split_node index_split breaks_split right_split node_party
#' @references \cite{Shannon, W. D., & Banks, D. (1999). Combining classification trees using MLE. Statistics in medicine, 18(6), 727-740.}
#' @references \cite{Chipman, H. A., George, E. I., & McCulloh, R. E. (1998). Making sense of a forest of trees. Computing Science and Statistics, 84-92.}
#' @references \cite{Sies, A. & Van Mechelen I. (Submitted). C443: An R-package to see a forest for the trees}
#' @examples 
#' require(MASS)
#' #Grow a forest of classification trees based on 10 bootstrap samples
#' forest <- growforest(Pima.tr, X=c("npreg", "glu", "bp", "skin", "bmi", "ped", "age"), 
#' Y ="type", ntrees = 10)
#' 
#' # Calculate similiarties between all pairs of trees in the forest
#' simmatrix <- similarities(fulldata = Pima.tr, treedata = forest[[2]], Y = rep("type", 10), 
#' X = c("npreg", "glu", "bp", "skin", "bmi", "ped", "age"), trees = forest[[1]], m = 1, weight = 0)
#' 
#' simmatrix2 <- similarities(fulldata = Pima.tr, treedata = forest[[2]], Y = rep("type",10), 
#' X = c("npreg", "glu", "bp", "skin", "bmi", "ped", "age"), trees = forest[[1]], m = 1, 
#' weight = 1, tol = c(3, 30, 10, 10, 5, 0.3, 10))



similarities<- function (fulldata, treedata, Y, X, trees,  m,  weight=NULL, A=NULL, tol=NULL, regime=FALSE ){
  if (m == 1){
    sim <- M1(fulldata, treedata, Y, X, trees, weight, A, tol)
  }
  if (m == 2){
    sim <- M2(fulldata, treedata, Y, trees, weight, A)
  }
  if (m == 3){
    sim <- M3(fulldata, treedata, Y, trees, A)
  }
  if (m == 4){
    sim <- M4(fulldata, treedata, Y, trees, A)
  }
  if (m == 5){
    sim <- M5(fulldata, treedata, Y, X, trees, A)
  }
  if(m==6){
    sim <- M6(fulldata, treedata, Y, X, trees, A, tol, regime)
  }
  return(sim)  
}


# Measure 1: number of splitvariables & splitpoints in common/total number of splitvariables largest tree
#If weights = 1, only splitvariables taken into account
M1 <- function (fulldata, treedata, Y, X, trees, weight,  A=NULL, tol=NULL){
  pamtrees <- as.pam(fulldata, treedata, Y, trees, A)[[1]]
  n <- length(pamtrees)
  simmatrix <- matrix(c(0), n, n)
  splits <- lapply(1:n, function(i) splitv(pamtrees[[i]]))
  s <- sapply(1:n, function (k) sapply(k:n, function (l) sim(splits[[k]], splits[[l]], X, weight, tol)))
  
  for (i in 1:n){
    si<- s[[i]]
    si[is.nan(si)] <- 1
    simmatrix[i, c(i:n)] <- si
  }
  
  ind <- lower.tri(simmatrix)   #make matrix symmetric
  simmatrix[ind] <- t(simmatrix)[ind] 
  return(simmatrix)
}

####### Measure 2: paths###########
M2<- function(fulldata, treedata, Y, trees, weight,  A=NULL){
  pamtrees <- as.pam(fulldata, treedata, Y, trees, A)[[1]]
  n <- length(pamtrees)
  simmatrix <- matrix(c(0), n, n)
  subs <- sapply (1:n, function (k) returnsubpaths(pamtrees[[k]])) #split up each path into all subpaths
  
  dis <- matrix(c(0), length(pamtrees), length(pamtrees))
  d <- sapply(1:n, function (s) sapply(s:n, function (j) dissim(subs[[s]], subs[[j]], weight) ))
  
  for (i in 1:n){
    dis[i, c(i:n)] <- d[[i]]
  }
  ind <- lower.tri(dis)   #make matrix symmetric
  dis[ind] <- t(dis)[ind] 
  sim <- 1 - dis
  return(sim)
}

####### Measure 3: classification labels agreement
M3<- function (fulldata, treedata, Y, trees, A=NULL){
  pamtrees <- as.pam(fulldata, treedata, Y, trees, A)
  g <- pamtrees[[2]]
  tree <- pamtrees[[1]]
  n <- length(tree)
  sim <- matrix(0, length(tree), length(tree))
  s <- sapply(1:n, function (s) sapply(s:n, function (j) mean(g[[s]]==g[[j]])))
  
  for (i in 1:n){
    sim[i,c(i:n)] <- s[[i]]
  }
  
  ind <- lower.tri(sim)   #make matrix symmetric
  sim[ind] <- t(sim)[ind] 
  return(sim)
}

###### Measure 4: Partition metric ################
M4<- function (fulldata, treedata, Y, trees,  A=NULL){
  pamtrees <- as.pam(fulldata, treedata, Y, trees, A)
  nodes <- pamtrees[[3]]
  n <- length(nodes)
  if( ! is.data.frame(treedata)){
    par<- lapply(1:n, function (s) part(treedata[[s]], nodes[[s]]))
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    clusterExport(cl, c("metr", "par", "treedata", "n"), envir=environment())
    si <- parSapply(cl, 1:n, function (s) sapply (s:n, function (j) metr(par[[s]], par[[j]], treedata[[s]])))
    stopCluster(cl)
  }else{
    par <- lapply(1:n, function (s) part(treedata,nodes[[s]]))
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    clusterExport(cl, c("metr", "par", "treedata", "n"), envir=environment())
    si <- parSapply(cl, 1:n, function (s) sapply (s:n, function (j) metr(par[[s]], par[[j]], treedata)))
    stopCluster(cl)
  }
  
  
  sim <- matrix(0, length(nodes), length(nodes))
  for (i in 1:n){
    sim[i, c(i:n)] <- si[[i]]
  }
  ind <- lower.tri(sim)   #make matrix symmetric
  sim[ind] <- t(sim)[ind] 
  return(sim)
}



########### measure 5: Disjunctive normal ########
M5 <- function(fulldata, treedata, Y, X, trees,  A=NULL){
  trees <- as.pam(fulldata, treedata, Y, trees,A)
  g <- trees[[2]]
  trees <- trees[[1]]
  n <-length(trees)
  s <- matrix(0, length(trees), length(trees))
  di <- lapply(1:length(trees), function (i) disjnorm(trees[[i]], g[[i]], fulldata, X))
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  clusterExport(cl, c("dis", "di", "n"), envir=environment())
  si <- parSapply(cl, 1:n, function (i) sapply (i:n, function (j) dis(di[[i]], di[[j]])))
  stopCluster(cl)
  
  for (i in 1:n){
    s[i, c(i:n)] <- si[[i]]
  }
  
  ind <- lower.tri(s)   #make matrix symmetric
  s[ind] <- t(s)[ind] 
  s<- round(s, digits=3)
  return(s)
}

M6 <- function (fulldata, treedata, Y, X, trees,  A=NULL, tol=NULL, regime=FALSE){
  pamtrees <- as.pam(fulldata, treedata, Y, trees, A)[[1]]
  n <- length(pamtrees)
  simmatrix <- matrix(c(0), n, n)
  
  s <- lapply(1:n, function (k) splitvsets(pamtrees[[k]][[2]]))
  # Voor elk blad met elk ander blad sim berekenen
  
  if(regime==TRUE){best<- lapply(1:n, function(k) which.max(c(mean(as.numeric(treedata[[k]][treedata[[k]][,A]==1,Y[k]])), mean(as.numeric(treedata[[k]][treedata[[k]][,A]==1,Y[k]])))))}
  if(regime==FALSE){best<- lapply(1:n, function(k) mean(as.numeric(treedata[[k]][,Y[k]]))> 0.5 )}
  sim <- sapply(1:n, function (k) sapply(k:n, function (l) simsets(s[[k]], s[[l]], X, tol, best[[k]],best[[l]])))
  
  for (i in 1:n){
    si<- sim[[i]]
    si[is.nan(si)] <- 1
    simmatrix[i, c(i:n)] <- si
  }
  
  ind <- lower.tri(simmatrix)   #make matrix symmetric
  simmatrix[ind] <- t(simmatrix)[ind] 
  return(simmatrix)
}



############################################################################################
####subfunctions###########################################################################
##########################################################################################
## y=vector of outcome variables, for each tree, data can be either list of data sets or 1 dataset
## ordata=original dataset
as.pam<- function(fulldata, treedata, Y, trees, A=NULL){
  g <- list(0)
  node <- list(0)
  n <- length(trees)
  for (i in 1:n){
    if(length(trees[[i]]) > 2){
      paths <- listrulesparty(x=trees[[i]])
      prednode <- predict(trees[[i]], newdata = fulldata, type = "node")
      predresp <- as.numeric(predict(trees[[i]], newdata =fulldata, type = "response"))-1
      frame <- matrix(c(0), length(unique(prednode)), 2)
      frame[, 1] <- sort(unique(prednode))
      frame[, 2] <- sapply(sort(unique(prednode)), function(k) mean(predresp[prednode == k]))   
      path1 <- paths[frame[, 2] == 1]
      trees[[i]] <- list(paths, path1)
      trees[[i]][[1]] <- sapply(1:length(trees[[i]][[1]]), function (k) strsplit(trees[[i]][[1]][k], " & "))
      trees[[i]][[2]] <- sapply(1:length(trees[[i]][[2]]), function (k) strsplit(trees[[i]][[2]][k], " & "))
      g[[i]] <- predresp
      node[[i]] <- as.numeric(prednode)
    }else{
      trees[[i]] <- list(numeric(0), numeric(0))
      if( ! is.data.frame(treedata)){
        y <- treedata[[i]][, Y[i]]
      } else{
        y <- treedata[, Y[i]]
      }
      if(is.null(A)){
        if(sum(y == 0, na.rm=T) > sum(y == 1, na.rm=T)){
          g1 <- 0
        } else{
          g1<- 1
        }
        g[[i]] <- rep(g1, length(y))
        node[[i]] <- rep(1, length(y))
      }
      
      if( ! is.null(A)){
        if( ! is.data.frame(treedata)){
          a <- treedata[[i]][, A]
        } else{
          a <- treedata[, A]
        }
        g[[i]] <- rep((which.max(c(mean(y[a == 0], na.rm=T), mean(y[a==1], na.rm=T))) -1),length(y))
        node[[i]] <- rep(1, length(y))
      }
    }
  }
  return(list(trees = trees, g = g, node = node))
}

### partykit:::.list.rules.party
listrulesparty <- function (x, i = NULL, ...) 
{
  if (is.null(i)) 
    i <- nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    ret <- sapply(i, listrulesparty, x = x)
    names(ret) <- if (is.character(i)) 
      i
    else names(x)[i]
    return(ret)
  }
  if (is.character(i) && !is.null(names(x))) 
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[, findx:ncol(dat), drop = FALSE]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) 
      dat <- x$data
  }
  else {
    fit <- NULL
    dat <- x$data
  }
  rule <- c()
  recFun <- function(node) {
    if (id_node(node) == i) 
      return(NULL)
    kid <- sapply(kids_node(node), id_node)
    whichkid <- max(which(kid <= i))
    split <- split_node(node)
    ivar <- varid_split(split)
    svar <- names(dat)[ivar]
    index <- index_split(split)
    
    if (is.factor(dat[, svar])) {
      if (is.null(index)) 
        index <- ((1:nlevels(dat[, svar])) > breaks_split(split)) + 
          1
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"", paste(slevels, 
                                               collapse = "\", \"", sep = ""), "\")", sep = "")
    }
    else {
      if (is.null(index)) 
        index <- 1:length(kid)
        breaks <- cbind(c(-Inf, breaks_split(split)), c(breaks_split(split), 
                                                      Inf))
      sbreak <- breaks[index == whichkid, ]
      right <- right_split(split)
      srule <- c()
      if (is.finite(sbreak[1])) 
        srule <- c(srule, paste(svar, ifelse(right, ">", 
                                             ">="), sbreak[1]))
      if (is.finite(sbreak[2])) 
        srule <- c(srule, paste(svar, ifelse(right, "<=", 
                                             "<"), sbreak[2]))
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(node_party(x))
  paste(rule, collapse = " & ")
}


#### Subfunctions M1 ######
splitv <- function (trees){
  paths <- trees[[1]]
  if(length(paths) > 0){
    pathsw <- paths
    leaves <- length(paths)
    l <- sapply(1:leaves, function (k) length(paths[[k]]))
    spaths <- lapply(1:leaves, function(k) sub(" <.", ':', paths[[k]]))
    spaths <- lapply(1:leaves, function(k) sub(" >..", ':', spaths[[k]]))
    paths <- lapply(1:leaves, function(k) sub(" <.*", '', paths[[k]]))
    paths <- lapply(1:leaves, function(k) sub(" >.*", '', paths[[k]]))
    
    splitvarsa <- unique(unlist(spaths))
    splitvars <- sub(":.*", '', splitvarsa)
    splitvarsa <- sub(".*:", '', splitvarsa)
    splitvarsa <- as.numeric(splitvarsa)
    
    nsplitvar1 <- length(splitvars)       # number of splits in tree i
  }else{
    splitvars <- NULL
    splitvarsa <- NULL
    nsplitvar1 <- 0
  }
  return(list(splitvars = splitvars, splitvarsa = splitvarsa, nsplitvar1 = nsplitvar1))
}



sim<- function (splits1, splits2, X, weight, tol){
  ### Only predictors
  if(weight == 0){
    common1 <- length(splits1[[1]][pmatch(splits1[[1]], splits2[[1]], nomatch = 0)]) #pmatch: no doubles
    total1 <- splits1[[3]] + splits2[[3]] - common1
    Jaccard <- common1 / total1
  }
  ### Also splitpoints
  if(weight == 1){
    # Create matrix That shows for each variable of the splits1 whether it is equal to each variable in splits2 
    same <- matrix(c(0), length(splits1[[1]]), length(splits2[[1]]))
    
    if(length(same > 0)){
      for (i in 1:length(splits1[[1]])){
        for(j in 1:length(splits2[[1]])){
          same[i, j] <- splits1[[1]][[i]] == splits2[[1]][[j]]
        }
      }
    
      # For those variables that are the same in both trees, put splitpoints in the matrix
      splitpoints <- matrix(c(rep(splits2[[2]], nrow(same))), nrow(same), ncol(same), byrow=T)
      s <- c(splitpoints)
      sm <- c(same)
      s[sm == 0] <- NA
      splitpoints <- matrix(s, nrow(same), ncol(same), byrow=F)
      
      # Get the right tolerance for each variables in splits1
      tsa <- splits1[[1]]
      t <- match(tsa,X)
      correct<- matrix(c(0), nrow(same), ncol(same))
      # Look whether splitpoint of splits2 is within tolerance zone of splits1
      for (i in 1:length(splits1[[1]])){
        d <- tol[t[i]]
        correct[i, ] <- findInterval(splitpoints[i, ], c(splits1[[2]][i] - d, splits1[[2]][i] + d) ) == 1
      }
      
      correct[is.na(correct)] <- 0
      g <- graph_from_incidence_matrix(correct)
      common<- max_bipartite_match(g)$matching_size
      
      total<- splits1[[3]] + splits2[[3]] - common
      Jaccard<- common / total  
    } else{
      if(length(splits1[[1]]) == 0 & length(splits2[[1]]) == 0){
        Jaccard <- 1
      } else{Jaccard <- 0}
    } 
  }
  return(Jaccard)
  
}



##### Subfuncties M2
returnsubpaths <- function(trees){
  paths <- trees[[1]]
  
  if(length(paths) > 0){
    leaves<- length(paths)
    l<- sapply(1:leaves, function (k) length(paths[[k]]))
    
    lastpaths <- lapply(1:leaves, function(k) sub("<.*", '', paths[[k]][l[k]])) # remove punctuation last attribut of path
    lastpaths <- lapply(1:leaves, function(k) sub(">.*", '', lastpaths[[k]]))
    for(j in 1:leaves){
      paths[[j]][l[j]]<- lastpaths[j]    
    }
    
    paths <- lapply(1:leaves, function(k) sub("<.*", '-', paths[[k]]))   #replace splitpoints with - or +
    paths <- lapply(1:leaves, function(k) sub(">.*", '+', paths[[k]]))
    upaths <- unique(paths)
    
    subpaths <- list(0)
    for(j in 1:length(upaths)){
      d<- list(0)
      d[[1]] <- upaths[[j]]
      if(length(upaths[[j]]) > 1){
        for (i in 2:length(upaths[[j]])){
          d[[i]] <- d[[i - 1]][- length(d[[i - 1]])]   #Split each path up into all subpaths
        }
      }
      subpaths[[j]] <- d
    }
    subpaths <- unlist(subpaths, recursive=FALSE)
    lastsubpaths <- lapply(1:length(subpaths), function(k) gsub("[[:punct:]]", '', subpaths[[k]][length(subpaths[[k]])]))  # remove punctuation from last attribute of each subpath
    for(j in 1:length(subpaths)){
      subpaths[[j]][length(subpaths[[j]])] <- lastsubpaths[j]    
    }
    
    subpaths<- unique(subpaths)
  }else{subpaths<- NULL}
  return(subpaths)
}


dissim<- function (subs1,subs2,weights){
  l <- sapply(1:length(subs1), function (k) length(subs1[[k]]))  # length of each subpath of each path
  d <- matrix(c(0))
  l2 <- sapply(1:length(subs2), function (k) length(subs2[[k]])) #length of each subpath of each path
  
  for (k in 1:(max(max(l), max(l2)))){
    a <- as.numeric(subs1[l == k] %in% subs2[l2 == k])  # number of subpaths of tree i in j
    b <- as.numeric(subs2[l2 == k] %in% subs1[l == k])  # number of subpaths of tree j in i
    
    if(weights==0){d[k] <- (length(a) + length(b)) - (sum(a) + sum(b))}   # if weighs 0 every dissimilar subpath weighted equally
    if(weights==1){d[k] <- 1 / k * ((length(a) + length(b)) - (sum(a) + sum(b)))} # if weights 1 every dissimilar subpath weighted by 1/length subpath
  }
  
  if(weights == 0){dis <- sum(d) / (length(l[l > 0]) + length(l2[l2 > 0]))}  # divide #dissimilar subpaths by maximum dissimilar subpaths
  wl <- sum(1 / l)
  wl2 <- sum(1 / l2)
  wl[is.infinite(wl)] <- 0
  wl2[is.infinite(wl2)] <- 0
  if(weights == 1){dis <- sum(d) / (wl+wl2)}  #divide weighted #dissimilar subpaths by maximum weighted # dissimilar subpaths
  
  dis[is.na(dis)] <- 0
  return(dis)
}


### Subfunctions M4
part <- function (data, tree1){
  t1 <- matrix(c(0), nrow(data), nrow(data))
  for(i in 1:nrow(data) - 1){
    t1[i,] <- tree1[i] == tree1  # for each observation with each other observation: Same leaf?
  }
  return(t1)
}

metr<- function (t1, t2, data){
  ind <- upper.tri(t1) 
  part<- sum(t1[ind] == t2[ind]) / choose(nrow(data), 2)
  return(part)
}

##### sUBFUNCTIONS M5
##############################################################
dis <- function (tree1,tree2){
  
  if( ! is.null(tree1) & ! is.null(tree2)){
   common <- matrix(c(0), length(tree1), length(tree2))
    
    for(i in 1:length(tree1)){
      for(j in 1: length(tree2)){
        common[i, j] <- length(unlist(tree1[i])[pmatch(unlist(tree1[i]), unlist(tree2[j]), nomatch = 0)])
      }
    }
   rows <- apply(common, 2, max)
   cols <- apply(common, 1, max)
   
   common <- min(sum(rows == length(tree1[[1]])), sum(cols == length(tree1[[1]]))) 
   total <- length(tree1) + length(tree2) - common
   s <- common / total
    
  }else if(is.null(tree1) & ! is.null(tree2) | ! is.null(tree1) & is.null(tree2)){
    s <- 0
  }else{
    s <- 1
  }
  
  return(s)
}



disjnorm <- function (trees, g, fulldata, X ){
  path <- trees[[2]]    #Paths that lead to leaf with predicted class 1
  
  if(length(path) > 0){
    leaves <- length(path)  # number of leaves
    l <- sapply(1:leaves, function (k) length(path[[k]]))     # for each path number of elements
    levels <- matrix(c(0), length(X), 2)
    
    for (i in 1:length(X)){
      levels[i, ] <- levels(fulldata[, X[i]])
    }
    
    levels <- cbind(X, levels)
    npaths <-path
    spaths <- lapply(1:leaves, function(k) sapply(1:l[k], function (j) sub(" %in%.*", '', path[[k]][j])))
   
    
    for(i in 1:length(X)){
      npaths <- lapply(1:leaves, function(k) sapply(1:l[k], function (j) sub(paste(levels[i,2]), '-', npaths[[k]][j])))
    }
    
    for(i in 1:length(X)){
      npaths <- lapply(1:leaves, function(k) sapply(1:l[k], function (j) sub(paste(levels[i,3]), '+', npaths[[k]][j])))
    }
    
    for(i in 1:length(X)){
       signs <- lapply(1:leaves, function (i) gsub("[^-+]+", "", npaths[[i]]))
    }
    
    paths <- lapply(1:leaves, function (i) paste(spaths[[i]], signs[[i]]))
    b <- sapply (1:length(spaths), function (s) X[pmatch(X, spaths[[s]], nomatch = 0) == 0])
    
    if(length(b) > 0){
      if(is.list(b)){
        for(i in 1:length(b)){
          now <- b[[i]]
          nowe <- now[now != ""]
          nowe <- sort(nowe)
          # Add those variables in plus and minus form to those paths
          if(length(nowe > 0)){
            a <- matrix(c(0), 2, length(nowe))
            a[1, ] <- c(paste(nowe, '+', sep=' '))
            a[2, ] <- c(paste(nowe, '-', sep=' '))
            grid <- do.call(expand.grid, split(a, col(a)))
            
            for (k in 1:nrow(grid)){
              paths <- c(paths, list(unlist(c(paths[i], paste(unlist(grid[k,c(1:ncol(grid))]))))))
            }
          } else {paths <- paths}
        }
      }else{
        now <- b
        nowe <- now[now != ""]
        nowe <- sort(nowe)
        # Add those variables in plus and minus form to those paths
        if(length(nowe > 0)){
          a <- matrix(c(0), 2, length(nowe))
          a[1, ] <- c(paste(nowe, '+', sep=' '))
          a[2, ] <- c(paste(nowe, '-', sep=' '))
          grid <- do.call(expand.grid, split(a, col(a)))
          paths <- lapply (1:nrow(grid), function (k) c(paths[[1]], paste(unlist(grid[k,c(1:ncol(grid))]))))
        } else {
          paths <- paths
          }
      }
      
    } else{paths <- paths}
    
    
    l <- sapply(1:length(paths), function (k) length(paths[[k]]))
    paths <- unique(paths)
    paths <- paths[lengths(paths) == max(l)]
  }else{
    if(g[1] == 0){
        paths <- NULL
    }else{
      a <- matrix(c(0), 2, length(X))
      a[1, ] <- c(paste(X, '+', sep=' '))
      a[2, ] <- c(paste(X, '-', sep=' '))
      grid <- do.call(expand.grid, split(a, col(a)))
      paths <- lapply (1:nrow(grid), function (k) paste(unlist(grid[k, c(1:ncol(grid))])))
    }
  }
  return(paths)
}


#### Subfunctions M6
splitvsets <- function (paths){
  splitvars<- list(length(paths))
  splitvarsa<- list(length(paths))
  nsplitvar1<- list(length(paths))
  
  #paths <- trees[[1]]
  if(length(paths) > 0){
    for(i in 1:length(paths)){
      pathsw <- paths[[i]]
      #leaves <- length(paths)
      #l <- sapply(1:leaves, function (k) length(paths[[k]]))
      #spaths <- sub(" <.", ':', unlist(pathsw, use.names=F))
     # spaths <- sub(" >..", ':', spaths)
      
     # splitvarsa1 <- unique(unlist(spaths))
      splitvarsa1<- unique(unlist(pathsw))
      splitvars1<- sub("\\-.*$", '', splitvarsa1)
      splitvars[[i]]<- sub("\\d.*$", '', splitvars1)
      
      splitvarsa1 <- sub(".*=", '', pathsw)
      splitvarsa1 <- sub(".*<", '', splitvarsa1)
      splitvarsa[[i]] <- c(as.numeric(splitvarsa1))
      
      #nsplitvar1[[i]] <- length(splitvars[[i]])       # number of splits in tree i
    }
  }else{
    splitvars <- NULL
    splitvarsa <- NULL
    #nsplitvar1 <- 0
  }
  return(list(splitvars = splitvars, splitvarsa = splitvarsa))
}



simsets<- function (paths1, paths2, X, tol,best1,best2){
  Jaccard<- matrix(c(0),length(paths1[[1]]),length(paths2[[1]]))
  if(length(paths1[[1]])&length(paths2[[1]])>0) {
    for (i in 1:length(paths1[[1]])){
      for (j in 1:length(paths2[[1]])){
        same <- matrix(c(0), length(paths1[[1]][[i]]), length(paths2[[1]][[j]]))
        
        if(length(same > 0)){
          for (k in 1:length(paths1[[1]][[i]])){
            for(m in 1:length(paths2[[1]][[j]])){
              same[k, m] <- paths1[[1]][[i]][[k]] == paths2[[1]][[j]][[m]]
            }
          }
          
          # For those variables that are the same in both trees, put splitpoints in the matrix
          splitpoints <- matrix(c(rep(paths2[[2]][[j]], nrow(same))), nrow(same), ncol(same), byrow=T)
          s <- c(splitpoints)
          sm <- c(same)
          s[sm == 0] <- NA
          splitpoints <- matrix(s, nrow(same), ncol(same), byrow=F)
          
          tsa<- paths1[[1]][[i]]
          tsa<- sub(" <.", '', tsa)
          tsa<- sub(" >..", '', tsa)
          t <- match(tsa,X)
          correct<- matrix(c(0), nrow(same), ncol(same))
          # Look whether splitpoint of splits2 is within tolerance zone of splits1
          for (h in 1:length(paths1[[1]][[i]])){
            d <- tol[t[h]]
            correct[h, ] <- findInterval(splitpoints[h, ], c(paths1[[2]][[i]][h] - d, paths1[[2]][[i]][h] + d) ) == 1
          }
          
          correct[is.na(correct)] <- 0
          g <- graph_from_incidence_matrix(correct)
          common<- max_bipartite_match(g)$matching_size
          
          total<- length(paths1[[1]][[i]]) + length(paths2[[2]][[j]]) - common
          Jaccard[i,j]<- common / total  
        }else{
          if(length(paths1[[1]][[i]]) == 0 & length(paths2[[1]][[j]]) == 0){
            Jaccard[i,j] <- 1
          } else{Jaccard[i,j] <- 0}
        }
      }
    }
    
    g <- graph_from_incidence_matrix(Jaccard, weighted=TRUE)
    common<- max_bipartite_match(g)$matching_weight
    sim<- common/min(length(paths1[[1]]), length(paths2[[2]]))
  } else{
    if(length(paths1[[1]]) == 0 & length(paths2[[1]]) == 0){
      if(best1==best2){
        sim <- 1 
      } else{sim<-0}
     
    } else{sim <- 0}
  }
  return(sim)
}

