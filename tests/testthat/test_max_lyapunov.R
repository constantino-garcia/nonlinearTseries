library(nonlinearTseries)
context("Maximal Lyapunov exponent")

test_that("estimates equal theoretical results", {
  
  # Lyapunov 
  h = henon(
    n.sample =  5000,
    n.transient = 100,
    a = 1.4,
    b = 0.3,
    start = c(0.63954883, 0.04772637),
    do.plot = FALSE
  )
  
  ts = h$x
  # message("\nCalculating maximal Lyapunov exponent of the Henon attractor\n")
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
  # message("expected: ",0.41," calculated: ",estimate(x),"\n")
  #relative error
  expect_true(abs((0.41 - estimate(x)) / 0.41) < 0.1)
  
  # Lorenz 
  l = lorenz(
    sigma = 16,
    rho = 40,
    beta = 4,
    start = c(-10, -11, 47),
    time =  seq(0, 100, by = 0.01),
    do.plot = FALSE
  )
  ts = l$x
  # message("\nCalculating maximal Lyapunov exponent of the Lorenz attractor(16,40,4)\n")
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    time.lag = 8,
    radius = 0.35,
    theiler.window = 30,
    min.neighs = 5,
    min.ref.points = length(ts),
    max.time.steps = 200,
    number.boxes = NULL,
    sampling.period = 0.01,
    do.plot = FALSE
  )
  # message("expected: ",1.37," calculated: ",estimate(x),"\n")
  #relative error
  expect_true(abs((1.37 - estimate(x)) / 1.37) < 0.1)
  
  # message("\nCalculating maximal Lyapunov exponent of the Lorenz attractor(10,45.92,4)\n")
  l = lorenz(
    sigma = 10,
    rho = 45.92,
    beta = 4,
    start = c(-10, -11, 47),
    time =  seq(0, 70, by = 0.01),
    do.plot = FALSE
  )
  ts = l$x
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
  # message("expected: ",1.50," calculated: ",estimate(x),"\n")
  #relative error
  expect_true(abs((1.50 - estimate(x)) / 1.50) < 0.1)
  
  # rossler
  r = rossler(
    a = 0.15,
    b = 0.2,
    w = 10,
    start = c(0, 0, 0),
    time = seq(0, 1000, 0.1),
    do.plot = FALSE
  )
  # message("\nCalculating maximal Lyapunov exponent of the Rossler attractor(0.15,0.2,10)\n")
  ts = r$x
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    time.lag = 12,
    radius = 0.1,
    theiler.window = 50,
    min.neighs = 5,
    min.ref.points = length(r),
    max.time.steps = 300,
    number.boxes = NULL,
    sampling.period = 0.1,
    do.plot = FALSE
  )
  # message("expected: ",0.090," calculated: ",estimate(x),"\n")
  # Rossler system oscillates and the estimation is difficult: use a range!
  expect_true((0.065 <= estimate(x)) && (estimate(x) <= 0.1))
  
  
  r = rossler(
    a = 0.2,
    b = 0.2,
    w = 5.7,
    start = c(0, 0, 0),
    time = seq(0, 1000, 0.1),
    do.plot = FALSE
  )
  ts = r$x
  # message("\nCalculating maximal Lyapunov exponent of the Rossler attractor(0.2,0.2,10)\n")
  x = maxLyapunov(
    time.series = ts,
    min.embedding.dim = 5,
    time.lag = 12,
    radius = 0.1,
    theiler.window = 50,
    min.neighs = 5,
    min.ref.points = length(r),
    max.time.steps = 300,
    number.boxes = NULL,
    sampling.period = 0.1,
    do.plot = FALSE
  )
  # message("expected: ",0.069," calculated: ",estimate(x),"\n")
  #relative error
  expect_true((0.065 <= estimate(x)) && (estimate(x) <= 0.1))
  
  
  # logistic map
  logmap = nonlinearTseries::logisticMap(do.plot = FALSE, start = 0.3)
  # message("\nCalculating maximal Lyapunov exponent of the logistic map(4)\n")
  x = maxLyapunov(
    time.series = logmap,
    min.embedding.dim = 2,
    time.lag = 1,
    radius = 0.01,
    theiler.window = 50,
    min.neighs = 5,
    min.ref.points = 5000,
    max.time.steps = 6,
    number.boxes = NULL,
    sampling.period = 1,
    do.plot = FALSE
  )
  # message("expected: ",0.693," calculated: ",estimate(x),"\n")
  #relative error
  expect_true(abs((0.693 - estimate(x)) / 0.693) < 0.1)
  
})
