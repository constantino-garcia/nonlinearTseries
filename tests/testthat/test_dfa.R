library(nonlinearTseries)
context("DFA")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(1)
  # generate/load data
  noise = rnorm(6000)
  brown_noise = function(len = 100) {
    cumsum(rnorm(len))
  }
  bnoise = brown_noise(6000)
  # load stored pink noise
  pnoise = readRDS("../testdata/pink_noise.RDS")
  # Compare estimates with theoretical results
  # Test I: white noise
  sm = 10
  sM = 1000
  dfa_result = dfa(
    time.series = noise,
    npoints = 10,
    window.size.range = c(sm, sM),
    do.plot = FALSE
  )
  # message("expected: 0.5--- estimate: ", estimate(dfa_result),"\n")
  expect_true(abs(0.50 - estimate(dfa_result)) / 0.50 < 0.1)
  
  # Test II: pink noise
  dfa_result = dfa(
    time.series = pnoise,
    npoints = 10,
    window.size.range = c(sm, sM),
    do.plot = FALSE
  )
  ## message("expected: 1--- estimate: ",dfa_result$alpha1,"\n")
  expect_true(abs(1 - estimate(dfa_result)) / 1 < 0.1)
  
  # Test III: brown noise
  dfa_result = dfa(
    time.series = bnoise,
    npoints = 10,
    window.size.range = c(sm, sM),
    do.plot = FALSE
  )
  expect_true(abs(1.50 - estimate(dfa_result)) / 1.50 < 0.1)
  ## message("expected: 1.5 --- estimate: ",dfa_result$alpha1,"\n")
  
  
})
