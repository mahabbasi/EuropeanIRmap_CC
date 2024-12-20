#' ------------------- Edit the trainDI to make it faster -------------
#' 
trainDI_editted <- function(model = NA,
                           train = NULL,
                           variables = "all",
                           weight = NA,
                           CVtest = NULL,
                           CVtrain = NULL,
                           method="L2",
                           useWeight = TRUE,
                           LPD = FALSE,
                           verbose = TRUE){
  
  # get parameters if they are not provided in function call-----
  if(is.null(train)){train = aoa_get_train(model)}
  if(length(variables) == 1){
    if(variables == "all"){
      variables = aoa_get_variables(variables, model, train)
    }
  }
  if(is.na(weight)[1]){
    if(useWeight){
      weight = aoa_get_weights(model, variables = variables)
    }else{
      message("variable are not weighted. see ?aoa")
      weight <- t(data.frame(rep(1,length(variables))))
      names(weight) <- variables
    }
  }else{
    
    
    weight <- user_weights(weight, variables)
    
  }
  
  # get CV folds from model or from parameters
  folds <-  aoa_get_folds(model,CVtrain,CVtest)
  CVtest <- folds[[2]]
  CVtrain <- folds[[1]]
  
  # check for input errors -----
  if(nrow(train)<=1){stop("at least two training points need to be specified")}
  
  # reduce train to specified variables
  train <- train[,na.omit(match(variables, names(train)))]
  
  train_backup <- train
  
  # convert categorial variables
  catupdate <- aoa_categorial_train(train, variables, weight)
  
  train <- catupdate$train
  weight <- catupdate$weight
  
  # scale train
  train <- scale(train)
  
  # make sure all variables have variance
  if (any(apply(train, 2, FUN=function(x){all(is.na(x))}))){
    stop("some variables in train seem to have no variance")
  }
  
  # save scale param for later
  scaleparam <- attributes(train)
  
  
  # multiply train data with variable weights (from variable importance)
  if(!inherits(weight, "error")&!is.null(unlist(weight))){
    train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})
  }
  
  
  # calculate average mean distance between training data
  
  # trainDist_avrg <- c()
  # trainDist_min <- c()
  
  if(method=="MD"){
    if(dim(train)[2] == 1){
      S <- matrix(stats::var(train), 1, 1)
    } else {
      S <- stats::cov(train)
    }
    S_inv <- MASS::ginv(S)
  }
  
  if (verbose) {
    message("Computing DI of training data...")
    pb <- txtProgressBar(min = 0,
                         max = nrow(train),
                         style = 3)
  }
  
  numCores <- 21 # Use one less core than available
  cl <-  makeCluster(numCores)
  registerDoParallel(cl)
  
  # Initialize variables
  trainDist_avrg <- vector("list", nrow(train))
  trainDist_min <- vector("list", nrow(train))
  
  # Parallel loop using foreach
  results <- foreach(i = seq_len(nrow(train)), .combine = rbind, .packages = c("CAST", 'FNN')) %dopar% {
    # Initialize local variables
    trainDist_avrg_local <- NA
    trainDist_min_local <- NA
    
    .alldistfun <- function(point, reference, method, sorted = TRUE,S_inv=NULL){
      
      if (method == "L2"){ # Euclidean Distance
        if(sorted){
          return(FNN::knnx.dist(reference, point, k = dim(reference)[1]))
        } else {
          return(FNN::knnx.dist(point,reference,k=1))
        }
      } else if (method == "MD"){ # Mahalanobis Distance
        if(sorted){
          return(t(sapply(1:dim(point)[1],
                          function(y) sort(sapply(1:dim(reference)[1],
                                                  function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) ))))))
        } else {
          return(t(sapply(1:dim(point)[1],
                          function(y) sapply(1:dim(reference)[1],
                                             function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) )))))
        }
      }
    }
    
    # Calculate distance to other training data
    trainDist <- matrix(.alldistfun(t(matrix(train[i, ])), train, method, sorted = FALSE, S_inv))
    trainDist[i] <- NA
    
    trainDist_avrg_local <- mean(trainDist, na.rm = TRUE)
    
    # Mask of any data that are not used for training for the respective data point (using CV)
    whichfold <- NA
    if (!is.null(CVtrain) & !is.null(CVtest)) {
      whichfold <- as.numeric(which(lapply(CVtest, function(x) { any(x == i) }) == TRUE))  # Index of the fold where i is held back
      if (length(whichfold) > 1) {
        stop("a datapoint is used for testing in more than one fold. currently this option is not implemented")
      }
      if (length(whichfold) != 0) {  # In case that a data point is never used for testing
        trainDist[!seq_len(nrow(train)) %in% CVtrain[[whichfold]]] <- NA  # Everything that is not in the training data for i is ignored
      }
      if (length(whichfold) == 0) {  # In case that a data point is never used for testing, the distances for that point are ignored
        trainDist <- NA
      }
    }
    
    if (length(whichfold) == 0) {
      trainDist_min_local <- NA
    } else {
      trainDist_min_local <- min(trainDist, na.rm = TRUE)
    }
    
    # Return results for this iteration
    c(trainDist_avrg_local, trainDist_min_local)
  }
  
  # Extract results from the combined matrix
  trainDist_avrg <- results[, 1]
  trainDist_min <- results[, 2]
  
  # Stop the parallel cluster
  stopCluster(cl)
  
  # If verbose is TRUE, display progress (this part can't be parallelized easily)
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = nrow(train), style = 3)
    for (i in seq_len(nrow(train))) {
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  
  # if (verbose) {
  #   close(pb)
  # }
  
  trainDist_avrgmean <- stats::median(trainDist_avrg,na.rm=TRUE)
  
  # Dissimilarity Index of training data -----
  TrainDI <- trainDist_min/trainDist_avrgmean
  
  
  # AOA Threshold ----
  # threshold_quantile <- stats::quantile(TrainDI, 0.75,na.rm=TRUE)
  # threshold_iqr <- (1.5 * stats::IQR(TrainDI,na.rm=T))
  # thres <- threshold_quantile + threshold_iqr
  thres <- mean(TrainDI) + 3*sd(TrainDI)
  # account for case that threshold_quantile + threshold_iqr is larger than maximum DI.
  if (thres>max(TrainDI,na.rm=T)){
    thres <- max(TrainDI,na.rm=T)
  }
  
  # note: previous versions of CAST derived the threshold this way:
  # thres <- grDevices::boxplot.stats(TrainDI)$stats[5]
  
  
  # calculate trainLPD and avrgLPD according to the CV folds
  if (LPD == TRUE) {
    if (verbose) {
      message("Computing LPD of training data...")
      pb <- txtProgressBar(min = 0,
                           max = nrow(train),
                           style = 3)
    }
    
    trainLPD <- c()
    for (j in  seq(nrow(train))) {
      
      # calculate  distance to other training data:
      trainDist      <- .alldistfun(t(matrix(train[j,])), train, method, sorted = FALSE, S_inv)
      DItrainDist <- trainDist/trainDist_avrgmean
      DItrainDist[j]   <- NA
      
      # mask of any data that are not used for training for the respective data point (using CV)
      whichfold <- NA
      if(!is.null(CVtrain)&!is.null(CVtest)){
        whichfold <- as.numeric(which(lapply(CVtest,function(x){any(x==j)})==TRUE)) # index of the fold where i is held back
        if(length(whichfold)>1){stop("a datapoint is used for testing in more than one fold. currently this option is not implemented")}
        if(length(whichfold)!=0){ # in case that a data point is never used for testing
          DItrainDist[!seq(nrow(train))%in%CVtrain[[whichfold]]] <- NA # everything that is not in the training data for i is ignored
        }
        if(length(whichfold)==0){#in case that a data point is never used for testing, the distances for that point are ignored
          DItrainDist <- NA
        }
      }
      
      #######################################
      
      if (length(whichfold)==0){
        trainLPD <- append(trainLPD, NA)
      } else {
        trainLPD <- append(trainLPD, sum(DItrainDist[,1] < thres, na.rm = TRUE))
      }
      if (verbose) {
        setTxtProgressBar(pb, j)
      }
    }
    
    if (verbose) {
      close(pb)
    }
    
    # Average LPD in trainData
    avrgLPD <- round(mean(trainLPD))
  }
  
  
  # Return: trainDI Object -------
  
  aoa_results = list(
    train = train_backup,
    weight = weight,
    variables = variables,
    catvars = catupdate$catvars,
    scaleparam = scaleparam,
    trainDist_avrg = trainDist_avrg,
    trainDist_avrgmean = trainDist_avrgmean,
    trainDI = TrainDI,
    threshold = thres,
    method = method
  )
  
  if (LPD == TRUE) {
    aoa_results$trainLPD <- trainLPD
    aoa_results$avrgLPD <- avrgLPD
  }
  
  class(aoa_results) = "trainDI"
  
  return(aoa_results)
}

################################################################################
# Helper functions
################################################################################
# Encode categorial variables

aoa_categorial_train <- function(train, variables, weight){
  
  # get all categorial variables
  catvars <- tryCatch(names(train)[which(sapply(train[,variables], class)%in%c("factor","character"))],
                      error=function(e) e)
  
  if (!inherits(catvars,"error")&length(catvars)>0){
    message("warning: predictors contain categorical variables. The integration is currently still under development. Please check results carefully!")
    
    for (catvar in catvars){
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
      train[,catvar]<-droplevels(train[,catvar])
      
      # then create dummy variables for the remaining levels in train:
      dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = train),
                           train)
      train <- data.frame(train,dvi_train)
      
      if(!inherits(weight, "error")){
        addweights <- data.frame(t(rep(weight[,which(names(weight)==catvar)],
                                       ncol(dvi_train))))
        names(addweights)<- colnames(dvi_train)
        weight <- data.frame(weight,addweights)
      }
    }
    if(!inherits(weight, "error")){
      weight <- weight[,-which(names(weight)%in%catvars)]
    }
    train <- train[,-which(names(train)%in%catvars)]
  }
  return(list(train = train, weight = weight, catvars = catvars))
  
  
}



# Get weights from train object


aoa_get_weights = function(model, variables){
  
  weight <- tryCatch(if(model$modelType=="Classification"){
    as.data.frame(t(apply(caret::varImp(model,scale=F)$importance,1,mean)))
  }else{
    as.data.frame(t(caret::varImp(model,scale=F)$importance[,"Overall"]))
  }, error=function(e) e)
  if(!inherits(weight, "error") & length(variables)>1){
    names(weight)<- rownames(caret::varImp(model,scale=F)$importance)
  }else{
    # set all weights to 1
    weight <- as.data.frame(t(rep(1, length(variables))))
    names(weight) = variables
    message("note: variables were not weighted either because no weights or model were given,
            no variable importance could be retrieved from the given model, or the model has a single feature.
            Check caret::varImp(model)")
  }
  
  #set negative weights to 0
  if(!inherits(weight, "error")){
    weight <- weight[,na.omit(match(variables, names(weight)))]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }
  return(weight)
  
  }



# check user weight input
# make sure this function outputs a data.frame with
# one row and columns named after the variables

user_weights = function(weight, variables){
  
  # list input support
  if(inherits(weight, "list")){
    # check if all list entries are in variables
    weight = as.data.frame(weight)
  }
  
  
  #check if manually given weights are correct. otherwise ignore (set to 1):
  if(nrow(weight)!=1  || !all(variables %in% names(weight))){
    message("variable weights are not correctly specified and will be ignored. See ?aoa")
    weight <- t(data.frame(rep(1,length(variables))))
    names(weight) <- variables
  }
  weight <- weight[,na.omit(match(variables, names(weight)))]
  if (any(weight<0)){
    weight[weight<0]<-0
    message("negative weights were set to 0")
  }
  
  return(weight)
  
}




# Get trainingdata from train object

aoa_get_train <- function(model){
  
  train <- as.data.frame(model$trainingData)
  return(train)
  
  
}


# Get folds from train object


aoa_get_folds <- function(model, CVtrain, CVtest){
  ### if folds are to be extracted from the model:
  if (!is.na(model)[1]){
    if(tolower(model$control$method)!="cv"){
      message("note: Either no model was given or no CV was used for model training. The DI threshold is therefore based on all training data")
    }else{
      CVtest <- model$control$indexOut
      CVtrain <- model$control$index
    }
  }
  ### if folds are specified manually:
  if(is.na(model)[1]){
    
    if(!is.null(CVtest)&!is.list(CVtest)){ # restructure input if CVtest only contains the fold ID
      tmp <- list()
      for (i in unique(CVtest)){
        tmp[[i]] <- which(CVtest==i)
      }
      CVtest <- tmp
    }
    
    if(is.null(CVtest)&is.null(CVtrain)){
      message("note: No model and no CV folds were given. The DI threshold is therefore based on all training data")
    }else{
      if(is.null(CVtest)){ # if CVtest is not given, then use the opposite of CVtrain
        CVtest <- lapply(CVtrain,function(x){which(!sort(unique(unlist(CVtrain)))%in%x)})
      }else{
        if(is.null(CVtrain)){ # if CVtrain is not given, then use the opposite of CVtest
          CVtrain <- lapply(CVtest,function(x){which(!sort(unique(unlist(CVtest)))%in%x)})
        }
      }
    }
  }
  return(list(CVtrain,CVtest))
}






# Get variables from train object

aoa_get_variables <- function(variables, model, train){
  
  if(length(variables) == 1){
    if(variables == "all"){
      if(!is.na(model)[1]){
        variables <- names(model$trainingData)[-which(names(model$trainingData)==".outcome")]
      }else{
        variables <- names(train)
      }
    }
  }
  return(variables)
  
  
}



.mindistfun <- function(point, reference, method, S_inv=NULL){
  
  if (method == "L2"){ # Euclidean Distance
    return(c(FNN::knnx.dist(reference, point, k = 1)))
  } else if (method == "MD"){ # Mahalanobis Distance
    return(sapply(1:dim(point)[1],
                  function(y) min(sapply(1:dim(reference)[1],
                                         function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) )))))
  }
}

.alldistfun <- function(point, reference, method, sorted = TRUE,S_inv=NULL){
  
  if (method == "L2"){ # Euclidean Distance
    if(sorted){
      return(FNN::knnx.dist(reference, point, k = dim(reference)[1]))
    } else {
      return(FNN::knnx.dist(point,reference,k=1))
    }
  } else if (method == "MD"){ # Mahalanobis Distance
    if(sorted){
      return(t(sapply(1:dim(point)[1],
                      function(y) sort(sapply(1:dim(reference)[1],
                                              function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) ))))))
    } else {
      return(t(sapply(1:dim(point)[1],
                      function(y) sapply(1:dim(reference)[1],
                                         function(x) sqrt( t(point[y,] - reference[x,]) %*% S_inv %*% (point[y,] - reference[x,]) )))))
    }
  }
}


