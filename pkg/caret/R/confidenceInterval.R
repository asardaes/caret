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
    
    if(min(B$t[ , 1L]) != max(B$t[ , 1L])) {
      metricCI <- boot::boot.ci(B, type = confType, L = L, conf = confLevel, ...)[[out_type]]
      
      if(confType != "norm") metricCI <- metricCI[-c(2:3)]
      metricCI <- c(B$t0, metricCI)
      
    } else metricCI <- c(B$t0, confLevel, B$t[1L], B$t[1L])
    
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
                                     number = 1000) {
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
        if(is.null(confGamma) || is.na(confGamma) || !is.numeric(confGamma)) 
          stop("Gamma was not estimated or was not provided. Please set its value manually.")
        L <- NULL
      }
      
      ## Final estimate given by training might not be simply the average value
      originalEstimate <- merge(object$results, object$bestTune)[[metric]]
      
      ## in case of trControl$returnResamp = "all"
      object <- merge(object$resample, object$bestTune)[[metric]]
      
      metricCI <- NextMethod("confidenceInterval",
                             confLevel = confLevel, confType = confType, confGamma = confGamma,
                             L = L, metric = metric, 
                             sampleSize = sampleSize, subsampleSizes = subsampleSizes, ...)
      
    } else stop("Resample results are not available.")
    
    metricCI[1L] <- originalEstimate
    metricCI
    
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
      
      B <- boot::boot(testdata, statistic = statistic, stype = "i", R = number, ...)
      
      metricCI <- boot::boot.ci(B, conf = confLevel, type = confType, ...)[[out_type]]
      if(confType != "norm") metricCI <- metricCI[-c(2:3)]
      metricCI <- c(B$t0, metricCI)
      
      names(metricCI) <- c(object$metric, "ConfLevel", "Lower", "Upper")
      attr(metricCI, "confType") <- confType
      
      metricCI
      
    } else {
      ## L confidence interval
      
      if(!is.numeric(confGamma) && !is.character(confGamma)) 
        stop("Please set gamma value manually or choose between 'range' and 'quantile' methods.")
      else if(is.character(confGamma))
        confGamma <- match.arg(confGamma, c("range", "quantile"))
      
      testdata <- data.frame(obs = newoutcome, pred = predict(object, newdata, ...))
      
      samp_size <- round(nrow(testdata) ^ (2/3))
      
      if(is.character(confGamma)) {
        sub_sizes <- round(nrow(testdata) ^ (1 / 2 * ((1 + (1:20) / (20 + 1)))))
        
        if(any(sub_sizes > samp_size)) {
          sub_sizes <- round(seq(from = samp_size, to = 1, length.out = 21L))
          sub_sizes <- sub_sizes[-21L]
        }
        
        lev <- levels(newoutcome)
        
        subsamples <- t(sapply(1:number, function(dummy) {
          ids <- lapply(sub_sizes, function(ss) sample(nrow(testdata), ss))
          sapply(ids, function(id) {
            subdf <- testdata[id, , drop = FALSE]
            object$control$summaryFunction(subdf, lev = lev, model = object)[metric]
          })
        }))
        
        if(confGamma == "range") {
          q1 <- t(apply(subsamples, 2, quantile, probs = seq(from = 0.25, to = 0.01, length.out = 10)))
          q2 <- t(apply(subsamples, 2, quantile, probs = seq(from = 0.75, to = 0.99, length.out = 10)))
          y_ij <- log(q2 - q1)
          
        } else {
          q0 <- t(apply(subsamples, 2, quantile, probs = seq(from = 0.75, to = 0.95, length.out = 15)))
          y_ij <- log(q0)
        }
        
        y_ij[is.infinite(y_ij)] <- NA
        
        y_i. <- rowMeans(y_ij, na.rm = TRUE)
        y_bar <- mean(y_ij, na.rm = TRUE)
        log_bin <- log(sub_sizes)
        log_bar <- mean(log_bin)
        confGamma <- sum((log_bin - log_bar)^2) # denominator
        confGamma <- (-sum((y_i. - y_bar) * (log_bin - log_bar))) / confGamma
      }
      
      object <- sapply(1:number, function(dummy) {
        id <- sample(nrow(testdata), samp_size)
        subdf <- testdata[id, , drop = FALSE]
        object$control$summaryFunction(subdf, lev = lev, model = object)[metric]
      })
      
      NextMethod("confidenceInterval",
                 confLevel = confLevel, confType = confType, confGamma = confGamma,
                 metric = metric, sampleSize = nrow(newdata), subsampleSizes = rep(samp_size, length(object)), ...)
    }
  }
}
