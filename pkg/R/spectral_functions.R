# @export
hamming = function(n, alpha = 0.54, beta = 0.46){
  alpha - beta * cospi (2 * (0:(n - 1)) / (n - 1))
}

# @export
hanning = function(n) {
  0.5 * (1 - cospi( 2* (0:(n - 1)) / (n - 1)))
}

# @export
blackman = function(n, alpha = 0.16) {
  a0 = (1 - alpha) / 2
  a1 = 1 / 2
  a2 = alpha / 2
  a0 - a1 * cospi(2 * (0:(n - 1)) / (n - 1)) + a2 * cospi(4 * (0:(n - 1)) / (n - 1)) 
}



#' @import zoo
# @export
spectrogram = function(x, ...){
  UseMethod("spectrogram")
}

# @export
spectrogram.default = function(x, width, by, pad = 0) {
  if (width < 2) {
    stop("width must be >= 2")
  }
  FUN = function(x, pad) { 
    n = width + pad
    final_len = ifelse(n %% 2 == 0, n / 2 + 1, (n + 1) / 2)
    (Mod(fft(c(x * hamming(width), rep(0,pad)))) ^ 2)[1:final_len]
  }
  spec = t(rollapply(x - mean(x), width = width, by = by,
                   FUN = FUN, pad = pad))
  time_spec = seq((width + 1) / 2, (width + 1) / 2 + (ncol(spec) - 1)* by, by = by)
  freq_spec = seq(0, 0.5, length = nrow(spec))
  # Create class
  class(spec) = c("spectrogram", class(spec))
  attr(spec,"freq") = freq_spec
  attr(spec, "time") = time_spec
  attr(spec, "width") = width
  attr(spec, "by") = by
  spec
}

# @export
spectrogram.ts = function(x, width, by, pad = 0) {
  spec = spectrogram(as.numeric(x), width = width, by = by, pad = pad)
  attr(spec,"freq") = attr(spec,"freq") * frequency(x)
  attr(spec, "time") = attr(spec, "time") / frequency(x) 
  spec
}

# @export
spectrogram.zoo = function(x, width, by, pad = 0) {
  spectrogram.ts(as.ts(x), width = width, by = by, pad = pad)
}


# @export
plot.spectrogram = function(x, color.palette = terrain.colors, 
                            xlab = "Time", ylab = "Frequency(Hz)",
                            main = "Spectrogram", 
                            zlim = range(x, finite = TRUE),
                            trim = FALSE, ...) {
  if (trim) {
    x[which(x < zlim[[1]], arr.ind = T)] = zlim[[1]]
    x[which(x > zlim[[2]], arr.ind = T)] = zlim[[2]]
  }
  filled.contour(attr(x, "time"), attr(x, "freq"), t(x),
                 color.palette = color.palette, 
                 xlab = xlab, ylab = ylab, main = main,
                 zlim = zlim, ...)
}

# @export
time.spectrogram = function(x) {
  attr(x, "time")  
}

# @export
freq = function(x){
  UseMethod("freq")
}

# @export
freq.spectrogram = function(x){
  attr(x, "freq")  
}
