library(nonlinearTseries)
context("Embedding dimension calculation")

test_that("estimates equal theoretical results", {
  skip_on_cran()
  set.seed(1)
  h = henon(do.plot = FALSE,
            start = c(0.4, 0.234),
            n.sample = 3000)
  l = lorenz(
    sigma = 10,
    rho = 28,
    beta = 8 / 3,
    start = c(-10, -11, 47),
    time =  seq(0, 30, by = 0.01),
    do.plot = FALSE
  )
  ik = ikedaMap(
    n.sample = 3000,
    n.transient = 100,
    start = c(0.5341, 0.278),
    do.plot = FALSE
  )
  
  #results expected: 4,2,4,none
  x = estimateEmbeddingDim(
    ik$x,
    time.lag = 1,
    max.embedding.dim = 6,
    threshold = 0.9,
    do.plot = FALSE
  )
  expect_equal(4, x)
  
  x = estimateEmbeddingDim(
    h$x,
    time.lag = 1,
    max.embedding.dim = 6,
    threshold = 0.9,
    do.plot = FALSE
  )
  expect_equal(2, x)
  
  x = estimateEmbeddingDim(
    l$x,
    time.lag = 15,
    max.embedding.dim = 6,
    threshold = 0.9,
    do.plot = FALSE
  )
  expect_equal(4, x)
  
})
