library(nonlinearTseries)
context("Correlation dimension port")

test_that("correlation dimension port is correct", {
  tolerance = 1e-2
  results = readRDS('../testdata/corrDimPort.RDS')
  # Henon I
  set.seed(1)
  h=henon(n.sample = 3000,n.transient = 100, a = 1.4, b = 0.3, 
          start = c(0.73954883, 0.04772637), do.plot = FALSE)
  
  ts = h$x
  mmin = 2
  mmax = 5
  time.lag = 1
  rmin = 10 ^ -2.35
  rmax = 10 ^ -2
  np = 100
  theiler.window = 5
  
  xr = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  expect_equal(results[[1]], xr, tolerance=tolerance)
  
  # Henon II (different parameters)
  rmin=0.001
  rmax=0.01
  h=henon(n.sample =  3000,n.transient = 100, a = 1.2, b = 0.3, 
          start = c(1,1), do.plot = FALSE)
  
  ts = h$x
  xr = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  expect_equal(results[[2]], xr, tolerance=tolerance)
  
  # Lorenz I
  l = lorenz(sigma = 10, rho = 28, beta = 8 / 3,
             start = c(-10, -11, 47),
             time =  seq(0, 70, by = 0.01),
             do.plot = FALSE
  )
  ts = l$x
  mmin = 3
  mmax = 6
  time.lag = 10
  rmin = 10 ^ -0.6
  rmax = 1
  np = 100
  theiler.window = 100
  xr = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  expect_equal(results[[3]], xr, tolerance=tolerance)
  
  
  # Rossler I 
  r = rossler(
    a = 0.15,
    b = 0.2,
    w = 10,
    start = c(0, 0, 0),
    time = seq(0, 2000, 0.1),
    do.plot = FALSE
  )
  ts = r$x
  mmin = 2
  mmax = 5
  time.lag = 12
  rmin = 10 ^ -0.4
  rmax = 10 ^ -0.2
  np = 100
  theiler.window = 100
  xr = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  expect_equal(results[[4]],xr,tolerance=tolerance)
  
  # Logistic map
  logmap = logisticMap(
    r = 3.5699456,
    n.sample = 5000,
    n.transient = 500,
    do.plot = FALSE
  )
  ts = logmap
  xr = corrDim(
    time.series = ts,
    min.embedding.dim = 2,
    max.embedding.dim = 4,
    time.lag = 1,
    min.radius = 10 ^ -5,
    max.radius = 10 ^ -3.5,
    n.points.radius = 10,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  expect_equal(results[[5]], xr, tolerance = tolerance)
  
})
