#' L Confidence Interval
#' 
#' Calculate a confidence interval based on subsampling.
#' 
#' The confidence interval, denoted here as LCI, is based on Politis and Romano (1994).
#' 
#' It requires almost no additional overhead for the calculations, but it is not entirely
#' nonparametric. The normalization exponent \code{confGamma} should be set according to the
#' expected sampling distribution of the statistic (in this case \code{metric}). A value of 0.5
#' would assume normality. It can be estimated during training using the range or quantile methods
#' discussed in Bertail, Politis and Romano (1999), which requires more memory since it must create
#' a matrix with 20 columns and (number of replications * number of tunning parameter combinations)
#' rows. In general, the quantile method is only valid for metrics that are strictly positive.
#' 
#' If the resulting confidence interval is \code{NULL}, it means no cross-validation was done, or it
#' is unsupported (e.g \code{oob}). If it has infinite values, it could be due to:
#' 
#' \itemize{
#'   \item Not enough replications were performed. Try increasing \code{number} or \code{repeats}
#'     appropriately in \code{\link{trainControl}}.
#'   \item The confidence level (\code{confLevel}) is too high. Could be solved with the above, or
#'     by decreasing the confidence level.
#'   \item The chosen metric is not very smooth. Look at its histogram using the \code{resample}
#'     element of \code{\link{train}}'s output.
#' }
#' 
#' If a custom \code{summaryFunction} is provided, please make sure its output has at least one
#' numeric element with the same name as the \code{metric} specified in \code{\link{train}},
#' otherwise it might not be detected correctly.
#' 
#' @note 
#' 
#' The resampling results of training need to be available for this calculations to be feasible.
#' These interval has asymptotic convergence, so the more replications, the better.
#' 
#' The value given to \code{\link{train}} in \code{trControl$confLevel} is the one that is used for
#' the calculation during training. This can be changed when calling this function.
#' 
#' @references 
#' 
#' Dimitris N. Politis and Joseph P. Romano. ``Large sample confidence regions based on subsamples
#' under minimal assumptions''. The Annals of Statistics, pages 2031-2050, 1994.
#' 
#' Bertail, Patrice, Dimitris N. Politis, and Joseph P. Romano. ``On subsampling estimators with
#' unknown rate of convergence''. Journal of the American Statistical Association 94.446 (1999):
#' 569-579.
#' 
#' @param trainResult An object returned by \code{\link{train}}.
#' @param confLevel The desired confidence level between 0 and 1.
#' @param ... Further arguments for \code{\link[stats]{predict}} when \code{newdata} is provided.
#' @param confGamma The normalization exponent. To estimate, use "range" or "quantile", otherwise
#' provide the numeric value. See details.
#' @param newdata New data to use to calculate the interval.
#' @param newoutcome The outcomes corresponding to \code{newdata}.
#' @param number Number of subsamples for the calculation of the interval. Only relevant if
#'   \code{newdata} is provided.
#'   
#' @author Alexis Sarda
#' @importFrom stats predict
#' @export
#'   
lci <- function(trainResult, confLevel = trainResult$control$confLevel, ..., 
                confGamma = if(is.null(newdata)) trainResult$control$confGamma else "range",
                newdata = NULL, newoutcome = NULL, number = 100L)
{
  if(is.null(confLevel)) return(NULL)
  if(xor(is.null(newdata), is.null(newoutcome))) stop("New data and new outcome have to be provided together")
  
  metric <- trainResult$metric
  subsampleSizes <- lengths(trainResult$control$indexOut) ## replaced when gamma is re-estimated
  
  if(is.null(newdata)) {
    ## re-compute with original data, maybe with a new confLevel
    if(!is.numeric(confGamma) || confGamma <= 0) stop("confGamma must be provided and must be positive in this context")
    
    sampleSize <- nrow(trainResult$trainingData)
    
    if(!is.null(trainResult$resample)) {
      ## final estimate given by training might not be simply the average value
      originalEstimate <- merge(trainResult$results, trainResult$bestTune)[[metric]]
      
      ## in case of trControl$returnResamp = "all"
      samples <- merge(trainResult$resample, trainResult$bestTune)[[metric]]
      
      metricCI <- lciHelper(samples, confLevel, confGamma, sampleSize, subsampleSizes, metric)
      metricCI[1L] <- originalEstimate
      
    } else stop("Resample results are not available")
    
  } else {
    ## in case the confidence interval is based on new, hopefully unseen data...
    if(is.character(confGamma))
      confGamma <- match.arg(confGamma, c("range", "quantile"))
    else if(!is.numeric(confGamma) || confGamma <= 0)
      stop("confGamma must be provided and must be positive if not estimated")
    
    testdata <- data.frame(obs = newoutcome, 
                           pred = stats::predict(trainResult, newdata = newdata, ...))
    
    sampleSize <- round(nrow(testdata) ^ (2/3))
    
    if(is.character(confGamma)) {
      subsampleSizes <- round(nrow(testdata) ^ (1 / 2 * ((1 + (1:20) / (20 + 1)))))
      
      if(any(subsampleSizes > sampleSize)) {
        subsampleSizes <- round(seq(from = sampleSize, to = 1, length.out = 21L))
        subsampleSizes <- subsampleSizes[-21L]
      }
      
      lev <- levels(newoutcome)
      
      subsamples <- t(sapply(1L:number, function(dummy) {
        ids <- lapply(subsampleSizes, function(ss) sample(nrow(testdata), ss))
        sapply(ids, function(id) {
          subdf <- testdata[id, , drop = FALSE]
          trainResult$control$summaryFunction(subdf, lev = lev, model = trainResult)[metric]
        })
      }))
      
      confGamma <- estimateGamma(subsamples, subsampleSizes, confGamma)
    }
    
    samples <- sapply(1L:number, function(dummy) {
      id <- sample(nrow(testdata), sampleSize)
      subdf <- testdata[id, , drop = FALSE]
      trainResult$control$summaryFunction(subdf, lev = lev, model = trainResult)[metric]
    })
    
    metricCI <- lciHelper(samples, confLevel, confGamma, sampleSize, subsampleSizes, metric)
  }
  
  metricCI
}

lciHelper <- function(samples, confLevel, confGamma, sampleSize, subsampleSizes, metric) {
  Tn <- mean(samples)
  Tao_n <- sampleSize ^ confGamma
  Tao_b <- subsampleSizes ^ confGamma
  
  obj_std <- Tao_b * (samples - Tn)
  x <- seq(from = min(obj_std), to = max(obj_std), length.out = 1000L)
  
  Ln <- sapply(x, function(x) mean(obj_std <= x))
  ## plot(x, Ln) to check
  
  alpha <- (1 - confLevel) / 2
  
  ## this could return Inf
  suppressWarnings(c_n <- c(min(x[Ln >= alpha]), min(x[Ln >= 1 - alpha])))
  
  metricCI <- Tn + c_n / Tao_n
  
  metricCI <- c(Tn, confLevel, metricCI)
  
  names(metricCI) <- c(metric, "ConfLevel", "Lower", "Upper")
  
  metricCI
}

estimateGamma <- function(samples, sizes, method) {
  if(method == "range") {
    q1 <- t(apply(samples, 2L, quantile, probs = seq(from = 0.25, to = 0.01, length.out = 10)))
    q2 <- t(apply(samples, 2L, quantile, probs = seq(from = 0.75, to = 0.99, length.out = 10)))
    y_ij <- log(q2 - q1)
    
  } else {
    q0 <- t(apply(samples, 2L, quantile, probs = seq(from = 0.75, to = 0.95, length.out = 15)))
    y_ij <- log(q0)
  }
  
  y_ij[is.infinite(y_ij)] <- NA
  y_i. <- rowMeans(y_ij, na.rm = TRUE)
  y_bar <- mean(y_ij, na.rm = TRUE)
  log_bin <- log(sizes)
  log_bar <- mean(log_bin)
  confGamma <- sum((log_bin - log_bar)^2) # denominator
  (-sum((y_i. - y_bar) * (log_bin - log_bar))) / confGamma
}
