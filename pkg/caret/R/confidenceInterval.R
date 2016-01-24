confidenceInterval <- function(object, conf.level = 0.95, conf.type = "bca", ...) {
  UseMethod("confidenceInterval")
}

confidenceInterval.default <- function(object, conf.level = 0.95, conf.type = "bca", L = NULL, metric = "Metric", ...) {
  conf.type <- match.arg(conf.type, c("norm", "basic", "perc", "bca"))
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
      
    } else metricCI <- c(B$t0, NA, NA, NA)
  } else metricCI <- c(B$t0, NA, NA, NA)
  
  names(metricCI) <- c(metric, "ConfLevel", "Lower", "Upper")
  
  attr(metricCI, "conf.type") <- conf.type
  
  metricCI
}

confidenceInterval.train <- function(object,
                                     conf.level = object$control$conf.level,
                                     conf.type = object$control$conf.type, ...) {
  force(conf.level)
  force(conf.type)
  
  if(!is.null(object$resample) && !is.null(object$control$conf.level)) {
    L <- merge(object$empInf, object$bestTune)
    L <- as.numeric(L[ , grepl("^\\.obs", colnames(L))])
    
    metric <- object$metric
    
    ## in case of trControl$returnResamp = "all"
    object <- merge(object$resample, object$bestTune)[[metric]]
    
    NextMethod("confidenceInterval",
               conf.level = conf.level, conf.type = conf.type, 
               L = L, metric = metric, ...)
    
  } else stop("Resample results are not available or confidence level was set to NULL")
}
