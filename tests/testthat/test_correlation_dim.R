library(nonlinearTseries)
context("Correlation dimension")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(1)
  # Henon: 1.25 +- 0.02 (Grassberger and Procaccia 1983) -----------
  set.seed(1)
  ts = henon(
    n.sample = 5000,
    n.transient = 100,
    start = c(0.73954883, 0.04772637),
    do.plot = FALSE
  )$x
  x = corrDim(
    time.series = ts,
    min.embedding.dim = 2,
    max.embedding.dim = 5,
    time.lag = 1,
    min.radius = 0.02,
    max.radius = 0.1,
    n.points.radius = 50,
    do.plot = FALSE,
    theiler.window = 30,
    number.boxes = 100
  )
  expect_equal(estimate(x), 1.25, tolerance = 2 * 10 ^ -2)
  
  # Lorenz I: 2.05 ± 0.01 (Grassberger and Procaccia 1983) -----------
  ts = lorenz(
    start = c(-10, -11, 47),
    time =  seq(0, 150, by = 0.01),
    do.plot = FALSE
  )$x
  
  x = corrDim(
    time.series = ts,
    min.embedding.dim = 4,
    max.embedding.dim = 7,
    time.lag = 10,
    min.radius = 0.5,
    max.radius = 1,
    n.points.radius = 50,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  expect_equal(estimate(x), 2.05, tolerance = 1 * 10 ^ -2)
  
  # Logistic map: 0.500 ± 0.005 (Grassberger and Procraccia 1983) -----------
  ts = logisticMap(
    r = 3.5699456,
    n.sample = 5000,
    n.transient = 500,
    do.plot = FALSE
  )
  x = corrDim(
    time.series = ts,
    min.embedding.dim = 2,
    max.embedding.dim = 4,
    time.lag = 1,
    min.radius = 1e-5,
    max.radius = 0.00031,
    n.points.radius = 10,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  expect_equal(0.500, estimate(x), tolerance = 5 * 10 ^ -3)
  
  # Rossler: 1.991 + 0.065 (http://sprott.physics.wisc.edu/chaos/comchaos.htm.)
  r = rossler(
    a = 0.2,
    b = 0.2,
    w = 5.7,
    start = c(1, 1, 1),
    time = seq(0, 300, 0.01),
    do.plot = FALSE
  )
  ts = r$x
  mmin = 2
  mmax = 5
  time.lag = 30
  rmin = 0.15
  rmax = 0.4
  np = 20
  theiler.window = 200
  x = corrDim(
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
  expect_equal(1.99, estimate(x), tolerance = 5 * 10 ^ -2)
})
