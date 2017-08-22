library(nonlinearTseries)
context("Poincare map Port")

test_that("test Poincare map port", {
  set.seed(0)
  
  timeSeries = list()
  for (i in 1:3) {
    timeSeries[[(i - 1) * 3 + 1]] = henon(do.plot = FALSE, start = c(0.4, 0.234),
                                          n.sample = 6000)$x
    timeSeries[[(i - 1) * 3 + 2]] = lorenz( sigma = 10, rho = 28, beta = 8 / 3,
                                            start = c(-10, -11, 47),
                                            time =  seq(0, 50, by = 0.01),
                                            do.plot = FALSE)$x
    timeSeries[[(i - 1) * 3 + 3]] = ikedaMap(n.sample = 6000, n.transient = 100,
                                             start = c(0.5341, 0.278), 
                                             do.plot = FALSE)$x
  } 
  results=c()
  for (i in seq_along(timeSeries)) {
    ts = timeSeries[[i]]
    ed = sample(2:10, 1)
    tl = sample(1:10, 1)
    radius = runif(1, 0.75, 1.25)
    nhv = runif(ed, -1, 1)
    hp = runif(ed, -1, 1)
    results[[i]] = poincareMap(ts, ed, tl, 
                               normal.hiperplane.vector = nhv,
                               hiperplane.point = hp)
  }
 
  oldResults = readRDS("../testdata/port_poincare_results.RDS")  
  for (j in seq_along(oldResults)) {
    expect_equal(oldResults[[j]], results[[j]])
  }
})
