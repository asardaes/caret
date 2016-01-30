confidenceInterval <- function(object, confLevel = 0.95, confType = "L", confGamma = NULL, ...) {
  UseMethod("confidenceInterval")
}

confidenceInterval.default <- function(object, confLevel = 0.95, confType = "L", confGamma = NULL, ...,
                                       L = NULL, sampleSize, subsampleSizes, metric = "Metric") {
  confType <- match.arg(confType, c("norm", "basic", "perc", "bca", "L"))
  
  if(confType != "L") {
    if(is.null(L)) stop("Empirical influence values need to be provided in L for this confidence interval")
    out_type <- switch(confType, norm = "normal", basic = "basic", perc = "percent", bca = "bca")
    
    B <- list(t0 = mean(object), 
              t = as.matrix(object), 
              R = NROW(object), 
              call = "")
    
    if(min(B$t[,1]) != max(B$t[,1])) {
      metricCI <- tryCatch(boot::boot.ci(B, type = confType, L = L, conf = confLevel, ...)[[out_type]],
                           warning = function(w) w,
                           error = function(e) e)
      
      if (!inherits(metricCI, "condition")) {
        if(confType != "norm") metricCI <- metricCI[-c(2:3)]
        metricCI <- c(B$t0, metricCI)
        
      } else metricCI <- c(B$t0, confLevel, -Inf, Inf)
    } else metricCI <- c(B$t0, confLevel, -Inf, Inf)
    
  } else {
    if(length(object) != length(subsampleSizes))
      stop("Length mismatch between 'object' and 'subsampleSizes'")
    
    Tn <- mean(object)
    Tao_n <- sampleSize ^ confGamma
    Tao_b <- subsampleSizes ^ confGamma
    
    obj_std <- Tao_b * (object - Tn)
    x <- seq(from = min(obj_std), to = max(obj_std), length.out = 1000L)
    
    Ln <- sapply(x, function(x) mean(obj_std <= x))
    # plot(x, Ln) to check
    
    alpha <- (1 - confLevel) / 2
    
    # this could return Inf
    suppressWarnings(c_n <- c(min(x[Ln >= alpha]), min(x[Ln >= 1 - alpha])))
    
    metricCI <- Tn + c_n / Tao_n
    
    metricCI <- c(Tn, confLevel, metricCI)
  }
  
  names(metricCI) <- c(metric, "ConfLevel", "Lower", "Upper")
  
  attr(metricCI, "confType") <- confType
  
  metricCI
}

confidenceInterval.train <- function(object,
                                     confLevel = object$control$confLevel,
                                     confType = object$control$confType,
                                     confGamma = object$control$confGamma,
                                     ..., 
                                     newdata = NULL,
                                     newoutcome = NULL,
                                     bootNum = 1000) {
  force(confLevel)
  force(confType)
  force(confGamma)
  
  if(is.null(confLevel)) return(NULL)
  
  if(xor(is.null(newdata), is.null(newoutcome))) stop("New data and new outcomes have to be provided together")
  
  metric <- object$metric
  
  if(is.null(newdata)) {
    sampleSize <- nrow(object$trainingData)
    subsampleSizes <- lengths(object$control$indexOut)
    
    if(!is.null(object$resample)) {
      if(confType != "L") {
        if(is.null(object$empInf)) stop("The empirical influence values are not present in the object")
        
        L <- merge(object$empInf, object$bestTune)
        L <- as.numeric(L[ , grepl("^\\.obs", colnames(L))])
        
      } else {
        if(is.null(confGamma) || !is.numeric(confGamma)) 
          stop("Gamma was not estimated or was not provided. Please set its value manually.")
        L <- NULL
      }
      
      ## in case of trControl$returnResamp = "all"
      object <- merge(object$resample, object$bestTune)[[metric]]
      
      NextMethod("confidenceInterval",
                 confLevel = confLevel, confType = confType, confGamma = confGamma,
                 L = L, metric = metric, sampleSize = sampleSize, subsampleSizes = subsampleSizes, ...)
      
    } else stop("Resample results are not available.")
    
  } else {
    ## in case the confidence interval is based on new, hopefully unseen data...
    
    confType <- match.arg(confType, c("norm", "basic", "perc", "bca", "L"))
    
    if(confType != "L") {
      ## bootstrap confidence intervals
      
      out_type <- switch(confType, norm = "normal", basic = "basic", perc = "percent", bca = "bca")
      
      statistic <- function(df, ids, ...) {
        lev <- levels(newoutcome)
        model <- object
        object$control$summaryFunction(df[ids, , drop = FALSE], lev = lev, model = model)[metric]
      }
      
      testdata <- data.frame(obs = newoutcome, pred = predict(object, newdata, ...))
      
      B <- boot::boot(testdata, statistic = statistic, stype = "i", R = bootNum, ...)
      
      metricCI <- boot::boot.ci(B, conf = confLevel, type = confType, ...)[[out_type]]
      if(confType != "norm") metricCI <- metricCI[-c(2:3)]
      metricCI <- c(B$t0, metricCI)
      
      names(metricCI) <- c(object$metric, "ConfLevel", "Lower", "Upper")
      attr(metricCI, "confType") <- confType
      
      metricCI
      
    } else {
      ## L confidence interval
      
      if(is.null(confGamma) || !is.numeric(confGamma)) 
        stop("Gamma was not estimated or was not provided. Please set its value manually.")
      
      testdata <- data.frame(obs = newoutcome, pred = predict(object, newdata, ...))
      
      sub_sizes <- round(nrow(newdata) ^ (1 / 2 * ((1 + (1:1000) / 1001))))
      
      lev <- levels(newoutcome)
      
      object <- sapply(sub_sizes, function(sub) {
        id <- sample(nrow(testdata), sub)
        subdf <- testdata[id, , drop = FALSE]
        object$control$summaryFunction(subdf, lev = lev, model = object)[metric]
      })
      
      NextMethod("confidenceInterval",
                 confLevel = confLevel, confType = confType, confGamma = confGamma,
                 metric = metric, sampleSize = nrow(newdata), subsampleSizes = sub_sizes, ...)
    }
  }
}
