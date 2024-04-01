# Spectral index
# @param regression.range in normalized frequency units
# @export
spectralIndex = function(x, regression.range = c(0, 0.1), do.plot = TRUE, 
                         spans = NULL, kernel = NULL, taper = 0.1,
                         pad = 0, demean = TRUE, detrend = FALSE,
                         na.action = na.fail) {
  regression.range = regression.range * frequency(x)
  x_pgram  = stats::spec.pgram(x, spans = spans, kernel = kernel,
                               taper = taper, pad = pad,
                               demean = demean, detrend = detrend,
                               na.action = na.action, plot = FALSE)
  fit = lm(log(x_pgram$spec) ~ log(x_pgram$freq), 
           subset = 
             (x_pgram$freq >= regression.range[[1]] & 
              x_pgram$freq <= regression.range[[2]])
  )
  spectral_index = spectralIndexObject(fit, x_pgram)
  if (do.plot) {
    plot(spectral_index)
  }
  spectral_index$lm_data = fit
  spectral_index
}

spectralIndexObject = function(fit, spec) {
  if (!inherits(fit, "lm")) {
    stop("An object of class lm is required to build a spectralIndex object")
  }
  if (!inherits(spec, "spec")) {
    stop("An object of class spec is required to build a spectralIndex object")
  }
  fit_sum = summary(fit)
  spectral_index_fit = fit_sum$coefficients[2,]
  # change sign of the spectral index estimation since our model is
  # log(spec) ~ -beta * log(freq) 
  spectral_index_fit[[1]] = -spectral_index_fit[[1]]
  spectral_index_fit[[3]] = -spectral_index_fit[[3]]
  
  spectral_index_fit = list(spectral_index_fit = spectral_index_fit,
                            PSD_fit = list(freq = exp(fit$model$`log(x_pgram$freq)`),
                                           spec = exp(fit$fitted.values)),
                            spec_data = spec)
  class(spectral_index_fit) = "spectralIndex"
  spectral_index_fit
}

#' @export
plot.spectralIndex = function(x, type = "l", log = "xy", 
                              ylab = "spectrum", xlab = "normalized frequency",
                              ...) {
  spec = x$spec_data
  plot(spec$freq, spec$spec, type = type, log = log, ylab = ylab, xlab = xlab, ...)  
  lines(x$PSD_fit$freq, x$PSD_fit$spec, col = 2, lty = 2, lwd = 4)
  position = ifelse(x$spectral_index_fit[[1]] > 0, "topright", "bottomright")
  legend(position, legend = paste(expression(beta), "=", 
                                  round(x$spectral_index_fit[[1]], 3)),
         bty = "n")
}


getHurst = function(x) {
  UseMethod("getHurst")
}


#' @export
getHurst.spectralIndex = function(x) {
  beta = x$spectral_index_fit[[1]]
  Hs = c(H_fgn = NA, H_fbm = NA)
  
  if (beta > 1) {
    Hs[[2]] = (beta - 1) / 2
  } else {
    Hs[[1]] = (beta + 1) / 2  
  }
  Hs
}

getAlpha = function(x) {
  UseMethod("getAlpha")  
}

#' @export 
getAlpha.spectralIndex = function(x) {
  (x$spectral_index_fit[[1]] + 1) / 2
}



# # @export
# localSpectralIndex = function(x, do.plot = TRUE,  spans = NULL, 
#                               kernel = NULL, taper = 0.1,
#                               pad = 0, demean = TRUE, detrend = FALSE,
#                               na.action = na.fail, span = 0.75){
#   x_pgram  = stats::spec.pgram(x, spans = spans, kernel = kernel,
#                                taper = taper, pad = pad,
#                                demean = demean, detrend = detrend,
#                                na.action = na.action, plot = FALSE)
#   data = data.frame(x = log(x_pgram$freq),
#                     y = log(x_pgram$spec))
#   fit = loess(y ~ x, data = data,
#               span = span, degree = 1)
#   # new_x is uniformly sampled in log space
#   new_x = seq(min(data$x),
#               max(data$x),
#               length.out = length(data$x))
#   fit_pr = predict(fit,newdata = data.frame(x = new_x))
#   local_si = - diff(fit_pr) / diff(new_x)
#   if (do.plot) { 
#     par(mfrow=c(2,1))
#     plot(x_pgram$freq, x_pgram$spec, log="xy", type="l")
#     lines(x_pgram$freq, exp(fitted(fit)) , col = 2, lwd = 2)
#     #plot(head(x_pgram$freq,-1), diff_data$y,type="l")
#     plot(head(x_pgram$freq,-1),
#          local_si , col = 2, lwd = 2,type="l")
#     par(mfrow=c(1,1))
#   }
#   list(freq = head(new_x, -1), local_si = local_si)
# }
