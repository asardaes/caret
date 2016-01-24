\name{confidenceInterval}
\alias{confidenceInterval}
\title{Compute confidence intervals for supported objects}
\description{
Compute confidence intervals for a given metric using the values obtained during training. All intervals except "L" are calculated using \code{\link[boot]{boot.ci}}, but note that not all of the intervals included there are supported, nor recommended. See Details.
}
\usage{
confidenceInterval(object, confLevel = 0.95, confType = "bca", confGamma = 0.5, ...)
}
\arguments{
  \item{object}{
  For now, only \code{train} objects are supported.
}
  \item{confLevel}{
  The confidence level for the interval. Must be between 0 and 1.
}
  \item{confType}{
  The type of confidence interval. See Details.
}
  \item{confGamma}{
  The normalization exponent for the "L" interval. See Details.
}
  \item{\dots}{
  Ignored for now.
}
}
\details{
For \code{confType}, the supported interval types from the \code{boot} package's \code{\link[boot]{boot.ci}} are "norm", "basic", "perc" and "bca". The "L" interval is a custom implementation.

The recommended intervals are "L" (L-CI), which is based on Politis and Romano (1994), or "bca", which is the bias-corrected and accelerated (BCa) interval discussed in DiCiccio and Efron (1996).

L-CI requires almost no additional overhead for the calculations, but it is not entirely nonparametric. The normalization exponent \code{confGamma} should be set according to the expected sampling distribution of the static (in this case \code{metric}). A value of 0.5 would assume normality. See the cited article for more information.

For the non-L confidence intervals, empirical influence values have to be obtained, which means there is a need to create a matrix with one column for every observation in the data, and one row for every row in \code{tuneGrid} (which can be as much as \code{tuneLength} ^ number of tunning parameters depending on what was provided to \code{\link{train}}). If parallel computation is used, that number has to be multiplied by the amount of parallel workers. This can be significant in terms of both time and memory if the dataset is large and/or many replications are made. If \code{confType} is "L", then the empirical influence computations are not performed.

If the resulting confidence interval is \code{NULL}, it means no cross-validation was done, or it is unsupported (e.g \code{oob}). If it has infinite values, it could be due to:

\itemize{
  \item Not enough replications were performed. Try increasing \code{number} or \code{repeats} appropriately in \code{\link{trainControl}}.
  \item The confidence level (\code{confLevel}) is too high. Could be solved with the above, or by decreasing the confidence level.
  \item The chosen metric is not very smooth. Look at its histogram using the \code{resample} element of \code{\link{train}}'s output.
}

If a custom \code{summaryFunction} is provided, please make sure its output has at least one numeric element with the same name as the \code{metric} specified in \code{\link{train}}, otherwise it might not be detected correctly.
}
\value{
A named numeric vector with 4 elements: the chosen metric's average value, the confidence level of the interval, and the upper and lower values of said interval.
}
\references{
Thomas J DiCiccio and Bradley Efron. Bootstrap confidence intervals. Statistical science, pages 189-212, 1996.

Bradley Efron. Nonparametric estimates of standard error: the jackknife, the bootstrap and other methods. Biometrika, 68(3):589-599, 1981.

Dimitris N Politis and Joseph P Romano. Large sample confidence regions based on subsamples under minimal assumptions. The Annals of Statistics, pages 2031-2050, 1994.
}
\author{
Alexis Sarda
}
\note{
The resampling results of training need to be available for this calculations to be feasible. These intervals have asymptotic convergence, so the more replications, the better.

The value given to \code{link{train}} in \code{trControl$confLevel} is the one that is used for the calculation during training. This can be changed when calling this function, as well as \code{confType}, so long as the necessary information is available in the \code{train} object.

Using \code{confType = "L"} during trainig avoids calculating the empirical influence values for the observations, which is one of the most time-consuming steps of the confidence interval calculation. However, this means that you will not be able to extract any other type of confidence interval post-training. The L-CI can always be extracted if the resampling results are available, just set a desired \code{confLevel}.
}
\seealso{
\code{\link[boot]{boot.ci}}
}
