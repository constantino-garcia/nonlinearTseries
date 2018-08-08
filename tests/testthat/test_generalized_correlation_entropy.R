library(nonlinearTseries)
context("Generalized correlation entropy")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(1)
  q = 2
  # Henon
  h = henon(
    n.sample = 5000,
    n.transient = 100,
    a = 1.4,
    b = 0.3,
    start = c(0.5702737, 0.3953898),
    do.plot = FALSE
  )
  ts = h$x
  mmin = 2
  mmax = 5
  time.lag = 1
  rmin = 10 ^ -2.6
  rmax = 0.01
  np = 20
  theiler.window = 5
  # message("\nComputing Renyi entropy for the Henon attractor (1,4,0.3)\n")
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = q,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  # relative error
  # message("expected: 1.26  estimate: ",estimate(x),"\n")
  expect_true(abs((1.26 - estimate(x)) / 1.26) < 0.1)
  
  # Henon II 
  rmin = 0.001
  rmax = 0.01
  h = henon(
    n.sample = 5000,
    n.transient = 100,
    a = 1.2,
    b = 0.3,
    start = c(0.7902737, 0.3953898),
    do.plot = FALSE
  )
  ts = h$x
  # message("\nComputing Renyi entropy for the Henon attractor (1,2,0.3)\n")
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = q,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  # message("expected: 1.200  estimate: ",estimate(x),"\n")
  # relative error
  expect_true(abs((1.200 - estimate(x)) / 1.200) < 0.1)
  
  # Lorenz
  l = lorenz(
    sigma = 10,
    rho = 28,
    beta = 8 / 3,
    start = c(-10, -11, 47),
    time = seq(0, 100, by = 0.01),
    do.plot = FALSE
  )
  ts = l$x
  mmin = 3
  mmax = 6
  time.lag = 10
  rmin = 0.01
  rmax = 1
  np = 5
  theiler.window = 100
  # message("\nComputing Renyi entropy for the Lorenz attractor (10,28,8/3)\n")
  x = corrDim(
    time.series = ts,
    min.embedding.dim = mmin,
    max.embedding.dim = mmax,
    corr.order = q,
    time.lag = time.lag,
    min.radius = rmin,
    max.radius = rmax,
    n.points.radius = np,
    do.plot = FALSE,
    theiler.window = theiler.window,
    number.boxes = 100
  )
  # message("expected: 2.05  estimate: ",estimate(x),"\n")
  # relative error
  expect_true(abs((2.05 - estimate(x)) / 2.05) < 0.1)
  
  # Logistic map
  logmap = logisticMap(
    r = 3.5699456,
    n.sample = 5500,
    n.transient = 500,
    start = 0.68,
    do.plot = FALSE
  )
  ts = logmap
  # message("\nComputing Renyi entropy for the logistic map (3.5699456)\n")
  x = corrDim(
    time.series = ts,
    min.embedding.dim = 2,
    max.embedding.dim = 4,
    corr.order = q,
    time.lag = 1,
    min.radius = 10 ^ -5,
    max.radius = 0.01,
    n.points.radius = 10,
    do.plot = FALSE,
    theiler.window = 100,
    number.boxes = 100
  )
  # message("expected: 0.538  estimate: ",estimate(x),"\n")
  # relative error
  expect_true(abs((0.538 - estimate(x)) / 0.538) < 0.1)
  
})
