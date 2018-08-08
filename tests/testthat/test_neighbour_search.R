library(nonlinearTseries)
context("Neighbour search")

test_that("the box assisted algorithm yields the same results as the brute-force approach", {
  skip_on_cran()
  set.seed(1)
  eps = 0.4
  for (k in c(2, 4, 5, 10)) {
    data = readRDS(paste0("../testdata/neigh_search_vectors_", k, ".RDS"))
    # solution computed with a brute force approach
    sol = readRDS(paste0("../testdata/neigh_search_sol_", k,".RDS"))
    
    neighs = findAllNeighbours(data, eps)
    neighs = lapply(neighs,
                function(x) {
                  if (length(x) > 0) {
                    sort(x)
                  } else{
                    x
                  }
                })
    expect_equal(sol, neighs)
  }
})
