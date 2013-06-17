test.surrogate.fft = function(){
 
  tol = 10^-7
  cat("generating surrogate data...\n")
  cat("checking power spectrum and mean...\n")
  for (i in 1:3){
    original = arima.sim(n = 1000+i, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)));
    surrogate.data = FFTsurrogate(original)
    original.spectrum = abs(fft(original))^2
    surrogate.spectrum = abs(fft(surrogate.data))^2
    checkEqualsNumeric( surrogate.spectrum, original.spectrum, tolerance = tol)
    checkEqualsNumeric( mean(surrogate.data), mean(original), tolerance = tol)
  }
  
  
}