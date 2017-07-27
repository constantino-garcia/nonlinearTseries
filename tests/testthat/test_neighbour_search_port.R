library(nonlinearTseries)
context("Neighbour Search port")

test_that("ported code yields same results", {
  ntakens = 500
  nrepeat = 3
  for (radius in c(10 ^ seq(-2, log10(2), len = 10), 20)) {
    for (i in seq_len(nrepeat)) {
      embeddingD = sample(2:10, 1)
      itak = sample(1:ntakens, 1)
      takens = matrix(rnorm(ntakens * embeddingD), nrow = ntakens)
      number.boxes = sample(c(4, 10, 50, 100, 1000), 1)
      
      expect_equal(
        oldNonlinearTseries::findAllNeighbours(takens, radius, number.boxes),
        rcppFindAllNeighbours(takens, radius, number.boxes)
       )
      expect_equal(
        rcppNeighbourSearch(takens, itak, radius, number.boxes),
        oldNonlinearTseries::neighbourSearch(takens, itak, radius, number.boxes)
      )
    } 
  }
})
