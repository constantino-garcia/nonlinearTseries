library(nonlinearTseries)
context("Surrogate data generation")

test_that("the surrogate data has the same linear parameters of the original linear process", {
  skip_on_cran() 
  set.seed(1)
  tol = 10 ^ -7
  # message("generating surrogate data...\n")
  # message("checking power spectrum and mean...\n")
  arma_pars = readRDS("../testdata/arma_pars.RDS")
  for (i in seq_len(nrow(arma_pars))) {
    original = arima.sim(n = 1000, 
                         list(ar = as.numeric(arma_pars[i, 1:5]), 
                              ma = as.numeric(arma_pars[i, 6:10]))
    )
    surrogate.data = drop(FFTsurrogate(original,n.samples = 1))
    original.spectrum = as.numeric(abs(fft(original)) ^ 2)
    surrogate.spectrum = as.numeric(abs(fft(surrogate.data)) ^ 2)
    expect_equal( surrogate.spectrum, original.spectrum, tolerance = tol)
    expect_equal( mean(surrogate.data), mean(original), tolerance = tol)
  }
})
