library(nonlinearTseries)
context("Information dimension")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(1)
  # Henon
  ts = henon(n.sample =  5000, n.transient = 100,
             a = 1.4, b = 0.3, start =  c(0.8253681, 0.6955566),
             do.plot = FALSE)$x
  leps = infDim(ts, min.embedding.dim = 2, time.lag = 1,
                min.fixed.mass = 0.01, max.fixed.mass = 0.1,
                number.fixed.mass.points = 100, radius = 0.0005, 
                increasing.radius.factor = sqrt(2), number.boxes = 100,
                number.reference.vectors = 1000, theiler.window = 50,
                kMax = 300, do.plot = FALSE)
  expect_equal(1.24, estimate(leps, do.plot = FALSE), tolerance = 5 * 10 ^ -2)

  # Sinai map
  ts = sinaiMap(a = 0.3, n.sample = 5000, start = c(0.23489, 0.8923),
                do.plot = FALSE)$x
  leps = infDim(ts, min.embedding.dim = 2, time.lag = 1, min.fixed.mass = 0.001,
                max.fixed.mass = 0.1, number.fixed.mass.points = 500,
                radius = 0.05, increasing.radius.factor = sqrt(2),
                number.boxes = 100, number.reference.vectors = 1000,
                theiler.window = 200, kMax = 500, do.plot = FALSE)
  expect_equal(1.68, estimate(leps, do.plot = FALSE), 
               tolerance = 5 * 10 ^ -2)
  
})
