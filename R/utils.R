#' @importFrom graphics abline layout legend lines par plot plot.new filled.contour
#' @importFrom stats IQR as.formula dist fft lm median rnorm runif sd ts time as.ts frequency na.fail
#' @importFrom grDevices terrain.colors
NULL

vectorizePar = function(par, N, default=1:N){
  if (is.null(par)) {
    par = default
  }
  if (length(par) < N) {
    par = rep(par, length.out = N)
  }
  par
}
                        

propagateTakensAttr = function(x, takens){
  attr(x, "id") = attr(takens,"id")
  attr(x, "time.lag") = attr(takens, "time.lag")
  attr(x, "embedding.dim") = attr(takens, "embedding.dim")
  x
}

# @export
standardize = function(x){
  sd.x = sd(x, na.rm = T)
  if (sd.x == 0) {
    warning("sd(x) == 0, returning x")
    x.st = x
  } else {
    x.st = (x - mean(x, na.rm = T)) / sd.x
  }
  x.st
}

# @export
normalize = function(x, min.val = 0, max.val = 1){
  if (min.val >= max.val) {
    stop("min.val should be < max.val")
  }
  xmin = min(x, na.rm = T)
  xmax = max(x, na.rm = T)
  if (xmax == xmin) {
    warning("xmax == xmin, returnin min.val")
    # in this way we preserve class(x)
    x.norm = 0 * x + xmin 
  } else {
    x.norm = (max.val - min.val) * (x - xmin) / (xmax - xmin) + min.val
  }
  x.norm
}
