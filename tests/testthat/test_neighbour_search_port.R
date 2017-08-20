library(nonlinearTseries)
context("Neighbour Search port")

test_that("ported code yields same results", {
  ntakens = 500
  nrepeat = 3
  set.seed(1)
  results = readRDS('../testdata/neighbourSearchPort.RDS')
  counter = 1
  for (radius in c(10 ^ seq(-2, log10(2), len = 10), 20)) {
    for (i in seq_len(nrepeat)) {
      embeddingD = sample(2:10, 1)
      itak = sample(1:ntakens, 1)
      takens = matrix(rnorm(ntakens * embeddingD), nrow = ntakens)
      number.boxes = sample(c(4, 10, 50, 100, 1000), 1)
      
      xxx = results[[counter]] 
      counter = counter + 1
      expect_equal(
        xxx,
        findAllNeighbours(takens, radius, number.boxes)
      )
      xxx = results[[counter]] 
      counter = counter + 1
      expect_equal(
        neighbourSearch(takens, itak, radius, number.boxes),
        xxx)
    } 
  }
})
