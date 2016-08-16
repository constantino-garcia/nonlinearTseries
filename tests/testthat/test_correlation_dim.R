library(nonlinearTseries)
context("Correlation dimension")

test_that("estimates equal theoretical results", {
  # Henon: 1.25 +- 0.02 (Grassberger and Procaccia 1983) -----------
  set.seed(1)
  ts = henon(n.sample = 5000,
            n.transient = 100,
            start = c(0.73954883, 0.04772637),
            do.plot = FALSE)$x
  x = rcppCorrDim(time.series = ts,
              min.embedding.dim = 2,
              max.embedding.dim = 5,
              time.lag = 1,
              min.radius = 0.02,
              max.radius = 0.1,
              n.points.radius = 50,
              do.plot = FALSE,
              theiler.window = 30,
              number.boxes = 100)
  expect_equal(estimate(x), 1.25, tolerance = 2 * 10 ^ -2)
  
  # Lorenz I: 2.05 ± 0.01 (Grassberger and Procaccia 1983) -----------
  ts = lorenz(start = c(-10, -11, 47),
             time =  seq(0, 150, by = 0.01),
             do.plot = FALSE)$x
  
  x = rcppCorrDim(time.series = ts,
              min.embedding.dim = 4,
              max.embedding.dim = 7,
              time.lag = 10,
              min.radius = 0.5,
              max.radius = 1,
              n.points.radius = 50,
              do.plot = FALSE,
              theiler.window = 100,
              number.boxes = 100)
  expect_equal(estimate(x), 2.05, tolerance = 1 * 10 ^ -2)
  
  # Logistic map: 0.500 ± 0.005 (Grassberger and Procraccia 1983) -----------
  ts = logisticMap(r = 3.5699456, n.sample = 5000, 
                   n.transient = 500, do.plot = FALSE )
  x = rcppCorrDim(time.series = ts,
              min.embedding.dim = 2,
              max.embedding.dim = 4,
              time.lag = 1,
              min.radius = 1e-5,
              max.radius = 0.00031,
              n.points.radius = 10,
              do.plot = FALSE,
              theiler.window = 100,
              number.boxes = 100)
  expect_equal(0.500, estimate(x), tolerance = 5 * 10 ^ -3)
  
})


