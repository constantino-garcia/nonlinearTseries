library(nonlinearTseries)
context("RQA")

test_that("test rqa", {
  set.seed(0)
  
  skip("RQA not implemented yet")
  
  timeSeries = list()
  for (i in 1:3) {
    timeSeries[[(i - 1) * 3 + 1]] = henon(do.plot = FALSE, start = c(0.4, 0.234),
                                          n.sample = 3000)$x
    timeSeries[[(i - 1) * 3 + 2]] = lorenz( sigma = 10, rho = 28, beta = 8 / 3,
                                            start = c(-10, -11, 47),
                                            time =  seq(0, 30, by = 0.01),
                                            do.plot = FALSE)$x
    timeSeries[[(i - 1) * 3 + 3]] = ikedaMap(n.sample = 3000, n.transient = 100,
                                             start = c(0.5341, 0.278), 
                                             do.plot = FALSE)$x
  } 
  timeSeries = lapply(timeSeries, function(x) {x / sd(x)})
  results = list()
  for (i in seq_along(timeSeries)) {
    ts = timeSeries[[i]]
    ed = sample(2:10, 1)
    tl = sample(1:10, 1)
    radius = runif(1, 1e-5, 1e-3)
    lmin = sample(2:5, 1)
    dtb = sample(2:5, 1)
    results[[i]] = rqa(
      time.series = ts,
      embedding.dim = ed,
      time.lag = tl,
      radius = radius,
      lmin = lmin,
      do.plot = FALSE,
      distanceToBorder = dtb
    )
  }
  oldResults = readRDS("../testdata/port_rqa_results.RDS")  
  for (j in names(oldResults)) {
    expect_equal(oldResults[[j]], results[[j]])
  }
})
