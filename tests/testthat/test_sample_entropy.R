library(nonlinearTseries)
context("Sample entropy")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(1)
  corr.order = 2
  # Henon
  h = henon(
    n.sample = 15000,
    n.transient = 100,
    a = 1.4,
    b = 0.3,
    start = c(0.78, 0.8165),
    do.plot = FALSE
  )
  ts = h$x
  mmin = 2
  mmax = 10
  time.lag = 1
  rmin = 10 ^ -2.6
  rmax = 0.01
  np = 20
  theiler.window = 5
  
  
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = corr.order,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  se = sampleEntropy(x, do.plot = FALSE)
  #see Grassberger-estimation of the kolmogorov entropy from a chaotic signal
  # message("Expected K2 = ",0.325," Estimated = ",mean(se$sample.entropy["9",]),"\n")
  # calculate relative error
  expect_true((mean(se$sample.entropy["9", ]) - 0.325) / 0.325 < 0.1)
  
  # Logistic maps
  mmin = 2
  mmax = 5
  rmin = 10 ^ -5
  rmax = 0.01
  time.lag = 1
  ts = logisticMap(
    r = 3.5,
    n.sample = 5500,
    n.transient = 500,
    do.plot = FALSE,
    start = 0.5
  )
  rmin = 10 ^ -5
  rmax = 0.01
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = corr.order,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = 10,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  se = sampleEntropy(x, do.plot = FALSE)
  #the expected value is 0, so we can not use a relative error
  expect_true(mean(se$sample.entropy["4", ]) < 0.01)
  
  ts = logisticMap(
    r = 3.6,
    n.sample = 5500,
    n.transient = 500,
    do.plot = FALSE,
    start = 0.955
  )
  rmin = 10 ^ -5
  rmax = 0.01
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = corr.order,
    time.lag = 1,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = 10,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  se = sampleEntropy(x, do.plot = FALSE)
  
  ## message("Expected K2 = ",0.204," Estimated = ",mean(se$sample.entropy[4,]),"\n")
  # calculate relative error
  expect_true((mean(se$sample.entropy["4", ]) - 0.204) / 0.204 < 0.1)
  
  rmin = 10 ^ -5
  rmax = 0.01
  ts = logisticMap(
    r = 3.8,
    n.sample = 5500,
    n.transient = 500,
    do.plot = FALSE,
    start = 0.5
  )
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = corr.order,
    time.lag = 1,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = 10,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  se = sampleEntropy(x, do.plot = FALSE)
  ## message("Expected K2 = ",0.442," Estimated = ",mean(se$sample.entropy[4,]),"\n")
  # calculate relative error
  expect_true((mean(se$sample.entropy["4", ]) - 0.442) / 0.442 < 0.1)
  
})
