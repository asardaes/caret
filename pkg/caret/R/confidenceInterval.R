confidenceInterval <- function(object, conf.level = 0.95, conf.type = "bca", conf.gamma = 0, ...) {
  UseMethod("confidenceInterval")
}

confidenceInterval.default <- function(object, conf.level = 0.95, conf.type = "bca", 
                                       L = NULL, conf.gamma = 0, subsample.sizes, 
                                       metric = "Metric", ...) {
  conf.type <- match.arg(conf.type, c("norm", "basic", "perc", "bca", "L"))
  
  if(conf.type != "L") {
    if(is.null(L)) stop("Empirical influence values need to be provided in L for this confidence interval")
    out.type <- switch(conf.type,
                       norm = "normal",
                       basic = "basic",
                       perc = "percent",
                       bca = "bca")
    
    B <- list(t0 = mean(object), 
              t = as.matrix(object), 
              R = NROW(object), 
              call = "")
    
    if(min(B$t[,1]) != max(B$t[,1])) {
      metricCI <- tryCatch(boot::boot.ci(B, type = conf.type, L = L, conf = conf.level, ...)[[out.type]],
                           warning = function(w) w,
                           error = function(e) e)
      
      if (!inherits(metricCI, "condition")) {
        if(conf.type != "norm") metricCI <- metricCI[-c(2:3)]
        metricCI <- c(B$t0, metricCI)
        
      } else metricCI <- c(B$t0, conf.level, -Inf, Inf)
    } else metricCI <- c(B$t0, conf.level, -Inf, Inf)
    
  } else {
    if(length(object) != length(subsample.sizes))
      stop("Length mismatch between 'object' and 'subsample.sizes'")
    
    Tn <- mean(object)
    Tao_n <- length(object) ^ conf.gamma
    Tao_b <- subsample.sizes ^ conf.gamma
    
    obj_std <- Tao_b * (object - Tn)
    x <- seq(from = min(obj_std), to = max(obj_std), length.out = 1000L)
    
    Ln <- sapply(x, function(x) mean(obj_std <= x))
    # plot(x, Ln) to check
    
    alpha <- (1 - conf.level) / 2
    
    # this could return Inf
    suppressWarnings(c_n <- c(min(x[Ln >= alpha]), min(x[Ln >= 1 - alpha])))
    
    metricCI <- Tn + c_n / Tao_n
    
    metricCI <- c(Tn, conf.level, metricCI)
  }
  
  names(metricCI) <- c(metric, "ConfLevel", "Lower", "Upper")
  
  attr(metricCI, "conf.type") <- conf.type
  
  metricCI
}

confidenceInterval.train <- function(object,
                                     conf.level = object$control$confLevel,
                                     conf.type = object$control$confType,
                                     conf.gamma = object$control$confGamma,
                                     ...) {
  force(conf.level)
  force(conf.type)
  subsample.sizes <- lengths(object$control$indexOut)
  
  if(!is.null(object$resample) && !is.null(object$control$confLevel)) {
    if(conf.type != "L" && object$control$confType == "L") {
      stop("confType must be set to something different than 'L' during training to be able to calculate this interval")
      
    } else if(conf.type != "L") {
      L <- merge(object$empInf, object$bestTune)
      L <- as.numeric(L[ , grepl("^\\.obs", colnames(L))])
      
    } else L <- NULL
    
    metric <- object$metric
    
    ## in case of trControl$returnResamp = "all"
    object <- merge(object$resample, object$bestTune)[[metric]]
    
    NextMethod("confidenceInterval",
               conf.level = conf.level, conf.type = conf.type, 
               L = L, metric = metric, subsample.sizes = subsample.sizes, ...)
    
  } else stop("Resample results are not available or confidence level was set to NULL")
}
