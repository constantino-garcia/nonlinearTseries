library(nonlinearTseries)
context("Maximal Lyapunov port")

test_that("checks port correctness", {
  set.seed(1)
  tolerance = 5e-2
  # Henon system
  ts = henon(n.sample =  5000, n.transient = 100, a = 1.4, b = 0.3,
             start = c(0.63954883, 0.04772637),
             do.plot = FALSE
  )$x
 x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 2,
    max.embedding.dim = 4,
    time.lag = 1,
    radius = 0.001,
    theiler.window = 4,
    min.neighs = 2,
    min.ref.points = 500,
    max.time.steps = 10,
    number.boxes = NULL,
    do.plot = FALSE
  )
  xr = rcppMaxLyapunov(
    time.series = ts,
    min.embedding.dim = 2,
    max.embedding.dim = 4,
    time.lag = 1,
    radius = 0.001,
    theiler.window = 4,
    min.neighs = 2,
    min.ref.points = 500,
    max.time.steps = 10,
    number.boxes = NULL,
    do.plot = FALSE
  )
  expect_equal(estimate(x), estimate(xr), tolerance = tolerance)
  expect_equal(attributes(x), attributes(xr))
  expect_equal(dim(x$s.function), dim(xr$s.function))
  
  
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
    min.embedding.dim = 5,
    max.embedding.dim = 6,
    time.lag = 8,
    radius = 0.30,
    theiler.window = 100,
    min.neighs = 5,
    min.ref.points = length(ts),
    max.time.steps = 300,
    number.boxes = NULL,
    sampling.period = 0.01,
    do.plot = FALSE
  )
  xr = rcppMaxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    max.embedding.dim = 6,
    time.lag = 8,
    radius = 0.30,
    theiler.window = 100,
    min.neighs = 5,
    min.ref.points = length(ts),
    max.time.steps = 300,
    number.boxes = NULL,
    sampling.period = 0.01,
    do.plot = FALSE
  )
  expect_equal(estimate(x, c(0.5, 2)), 
               estimate(xr, c(0.5, 2)), tolerance = tolerance)
  expect_equal(attributes(x), attributes(xr))
  expect_equal(dim(x$s.function), dim(xr$s.function))
  
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
    max.embedding.dim = 10,
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
  xr = rcppMaxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    max.embedding.dim = 10,
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
  expect_equal(estimate(x), estimate(xr), tolerance = tolerance)
  expect_equal(attributes(x), attributes(xr))
  expect_equal(dim(x$s.function), dim(xr$s.function))
  
  
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
  xr = rcppMaxLyapunov(
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
  expect_equal(estimate(x), estimate(xr), tolerance = tolerance)
  expect_equal(attributes(x), attributes(xr))
  expect_equal(dim(x$s.function), dim(xr$s.function))
  
  # logistic map
  logmap = nonlinearTseries::logisticMap(do.plot = FALSE, n.sample = 10000, start = 0.3)
 x = maxLyapunov(
    time.series = logmap,
    min.embedding.dim = 2,
    max.embedding.dim = 7,
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
  xr = rcppMaxLyapunov(
    time.series = logmap,
    min.embedding.dim = 2,
    max.embedding.dim = 7,
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
  expect_equal(estimate(x), estimate(xr), tolerance = tolerance)
  expect_equal(attributes(x), attributes(xr))
  expect_equal(dim(x$s.function), dim(xr$s.function))
  
})
