library(nonlinearTseries)
context("Neighbour search")

# Regression test for issue #8 (Calls to Rf_error):
# nn.search uses the ANN library internally, which called Rf_error() on invalid
# inputs. This was replaced with Rcpp::stop() as required by RcppCore/Rcpp#1247.
# This test checks that errors from the ANN library are properly propagated as
# R conditions.
test_that("nn.search propagates ANN errors as R conditions", {
  set.seed(1)
  data <- matrix(rnorm(10), ncol = 2)  # 5 points
  ND <- nrow(data)
  NQ <- nrow(data)
  k <- 10L  # more neighbours than data points
  # Call the C++ wrapper directly, bypassing the R-level guard in nn.search,
  # to trigger the error inside the ANN library (kd_search.cpp).
  # It should surface as a catchable R error.
  expect_error(
    .Call("_nonlinearTseries_get_NN_2Set_wrapper",
          data, data, 2L, ND, NQ, k, 0.0, 1L, 0L, 1.0,
          integer(k * NQ), double(k * NQ)),
    regexp = "RANN"
  )
})

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
