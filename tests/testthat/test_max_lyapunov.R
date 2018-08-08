library(nonlinearTseries)
context("Maximal Lyapunov exponent")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(10)
  tolerance = 5e-2
  # Henon system
  ts = henon(n.sample =  5000, n.transient = 100, a = 1.4, b = 0.3,
             start = c(0.63954883, 0.04772637),
             do.plot = FALSE
  )$x
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 2,
    time.lag = 1,
    radius = 0.001,
    theiler.window = 4,
    min.neighs = 2,
    min.ref.points = 500,
    max.time.steps = 10,
    number.boxes = NULL,
    do.plot = FALSE
  )
  expect_equal(estimate(x), 0.41, tolerance = tolerance)
  
  # Lorenz 
  ts = lorenz(
    sigma = 16,
    rho = 40,
    beta = 4,
    start = c(-10, -11, 47),
    time =  seq(0, 100, by = 0.01),
    do.plot = FALSE
  )$x
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 4,
    max.embedding.dim = 6,
    time.lag = 8,
    radius = 0.40,
    theiler.window = 200,
    min.neighs = 5,
    min.ref.points = length(ts),
    max.time.steps = 400,
    number.boxes = NULL,
    sampling.period = 0.01,
    do.plot = FALSE
  )
  expected = 1.37
  expect_equal(estimate(x, c(0, 2.7), FALSE), expected, tolerance = tolerance, 
               scale = expected)
  
  # lorenz
  ts = lorenz(
    sigma = 10,
    rho = 45.92,
    beta = 4,
    start = c(-10, -11, 47),
    time =  seq(0, 70, by = 0.01),
    do.plot = FALSE
  )$x
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    time.lag = 8,
    radius = 0.5,
    theiler.window = 30,
    min.neighs = 5,
    min.ref.points = length(ts),
    max.time.steps = 200,
    number.boxes = NULL,
    sampling.period = 0.01,
    do.plot = FALSE
  )
  expect_equal(estimate(x), 1.50, tolerance = tolerance)
  
  
  # rossler
  ts = rossler(
    a = 0.15,
    b = 0.2,
    w = 10,
    start = c(0, 0, 0),
    time = seq(0, 3000, 0.1),
    do.plot = FALSE
  )$x
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    time.lag = 12,
    radius = 0.05,
    theiler.window = 50,
    min.neighs = 10,
    min.ref.points = 1000,
    max.time.steps = 300,
    number.boxes = NULL,
    sampling.period = 0.1,
    do.plot = FALSE
  )
  expect_equal(estimate(x), 0.065, tolerance = tolerance)
  
  # logistic map
  logmap = nonlinearTseries::logisticMap(do.plot = FALSE, n.sample = 10000, start = 0.3)
  x = maxLyapunov(
    time.series = logmap,
    min.embedding.dim = 2,
    time.lag = 1,
    radius = 0.01,
    theiler.window = 100,
    min.neighs = 10,
    min.ref.points = 10000,
    max.time.steps = 5,
    number.boxes = NULL,
    sampling.period = 1,
    do.plot = FALSE
  )
  expect_equal(estimate(x), 0.693, tolerance = tolerance)
  
})
