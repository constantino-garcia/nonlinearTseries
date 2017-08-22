library(nonlinearTseries)
context("nonlinear noise reduction port")

test_that("ported code yields same results", {
  set.seed(1)
  results = readRDS('../testdata/noiseReductionPort.RDS')
  counter = 1
  nrepeat = 3
  for (radius in c(10 ^ seq(-2, log10(2), len = 10), 20)) {
    for (i in seq_len(nrepeat)) {
      embeddingD = sample(2:10, 1)
      number.boxes = sample(c(4, 10, 50, 100, 1000), 1)
      ts = arima.sim(model = list(ar = runif(1, -1, 1),
                                  ma = runif(1, -1 , 1)), 
                     1000)
      cts1 = cts2 = ts
      expect_equal(
        results[[counter]],
        nonLinearNoiseReduction(cts2, embeddingD, radius)
      )
      counter = counter + 1
      # nonLinearNoiseReduction modifies directly the time
      # series in the C++ code, check that it doesn't affect the
      # R object
      expect_equal(cts2, ts)
    } 
  }
})
