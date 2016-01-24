confidenceInterval <- function(object, confLevel = 0.95, confType = "bca", confGamma = 0.5, ...) {
  UseMethod("confidenceInterval")
}

confidenceInterval.default <- function(object, confLevel = 0.95, confType = "bca", confGamma = 0.5, ...,
                                       L = NULL, subsampleSizes, metric = "Metric") {
  confType <- match.arg(confType, c("norm", "basic", "perc", "bca", "L"))
  
  if(confType != "L") {
    if(is.null(L)) stop("Empirical influence values need to be provided in L for this confidence interval")
    out.type <- switch(confType,
                       norm = "normal",
                       basic = "basic",
                       perc = "percent",
                       bca = "bca")
    
    B <- list(t0 = mean(object), 
              t = as.matrix(object), 
              R = NROW(object), 
              call = "")
    
    if(min(B$t[,1]) != max(B$t[,1])) {
      metricCI <- tryCatch(boot::boot.ci(B, type = confType, L = L, conf = confLevel, ...)[[out.type]],
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
    Tao_n <- length(object) ^ confGamma
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
                                     ...) {
  force(confLevel)
  force(confType)
  subsampleSizes <- lengths(object$control$indexOut)
  
  if(is.null(confLevel)) return(NULL)
  
  if(!is.null(object$resample)) {
    if(confType != "L") {
      if(is.null(object$empInf)) stop("The empirical influence values are not present in the object")
      
      L <- merge(object$empInf, object$bestTune)
      L <- as.numeric(L[ , grepl("^\\.obs", colnames(L))])
      
    } else L <- NULL
    
    metric <- object$metric
    
    ## in case of trControl$returnResamp = "all"
    object <- merge(object$resample, object$bestTune)[[metric]]
    
    NextMethod("confidenceInterval",
               confLevel = confLevel, confType = confType, 
               L = L, metric = metric, subsampleSizes = subsampleSizes, ...)
    
  } else stop("Resample results are not available.")
}
