library(nonlinearTseries)
context("Information dimension port")

test_that("ported code yields same results", {
  set.seed(1)
  results = readRDS('../testdata/infDimPort.RDS')
  ntakens = 500
  nrepeat = 2
  delta.embedding = 0:1
  counter = 1
  for (radius in c(10 ^ seq(-2, 0, len = 3))) {
    for (i in seq_len(nrepeat)) {
      min.embedding.dim = sample(2:4, 1)
      max.embedding.dim = min.embedding.dim + delta.embedding[[i]]
      min.fixed.mass = sample(c(1e-5,1e-4,1e-3), 1)
      max.fixed.mass = sample(c(5e-3,1e-2,5e-2), 1)
      increasing.radius.factor = 1 + sample(c(0.1,0.4, 0.6), 1)
      number.reference.vectors = sample(c(100,200,500), 1)
      kMax = sample(c(100,200,300), 1)
      time.lag = sample(1:10, 1)
      number.boxes = sample(c(4, 10, 50, 100), 1)
      ts = rnorm(2000)
      
      expect_equal(
        results[[counter]],
        rcppInfDim(ts, 
                   min.embedding.dim = min.embedding.dim,
                   max.embedding.dim = max.embedding.dim,
                   time.lag = time.lag,
                   min.fixed.mass = min.fixed.mass, max.fixed.mass = max.fixed.mass,
                   number.fixed.mass.points = 100, radius = radius, 
                   increasing.radius.factor = increasing.radius.factor,
                   number.boxes = 100,
                   number.reference.vectors = number.reference.vectors,
                   theiler.window = 50,
                   kMax = kMax, do.plot = FALSE)
        )
      counter = counter + 1
    }
  }
})
