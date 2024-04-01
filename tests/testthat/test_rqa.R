library(nonlinearTseries)
context("RQA")

test_that("empty RQA", {
  ts = seq(0, 10)
  ts_rqa = rqa(time.series = ts, embedding.dim = 2, time.lag = 1, radius = 0.1)
  rmatrix = as.matrix(ts_rqa$recurrence.matrix) 
  expect_true(all(rmatrix == diag(ncol(rmatrix))))
  expect_true(
    all(head(ts_rqa$diagonalHistogram, -1) == 0)
  )
})

test_that("full RQA", {
  ts = seq(0, 10)
  ts_rqa = rqa(time.series = ts, embedding.dim = 2, time.lag = 1, radius = 10, 
               lmin = 1)
  rmatrix = as.matrix(ts_rqa$recurrence.matrix) 
  expect_true(all(rmatrix == 1))
  # expected_diag = c(2,2,2,..., 1)
  expected_diag = rep(2, ncol(rmatrix))
  expected_diag[length(expected_diag)] = 1
  expect_equal(ts_rqa$diagonalHistogram, expected_diag)
  expected_vertical = rep(0, ncol(rmatrix))
  expected_vertical[ncol(rmatrix)] = ncol(rmatrix)
  expect_equal(ts_rqa$verticalHistogram, expected_vertical)
})


test_that("handmade RQA", {
  skip_on_cran()
  set.seed(1)
  # Test 1
  ts = readRDS("../testdata/rqa_ts_1.RDS")  
  handmade_result = readRDS("../testdata/rqa_results_1.RDS")
  rqa_result = rqa(
    time.series = ts,
    embedding.dim = 2,
    time.lag = 1,
    radius = 1.2,
    lmin = 2,
    do.plot = FALSE,
    distanceToBorder = 2
  )
  
  for (i in names(handmade_result)) {
    expect_equal(handmade_result[[i]], rqa_result[[i]])
  }
  
  # Test 2
  # Same ts as in Test 1
  handmade_result = readRDS("../testdata/rqa_results_2.RDS")
  rqa_result = rqa(
    time.series = ts,
    embedding.dim = 2,
    time.lag = 1,
    radius = 1.2,
    lmin = 3,
    vmin = 4,
    do.plot = FALSE,
    distanceToBorder = 3
  )
  
  for (i in names(handmade_result)) {
    expect_equal(handmade_result[[i]], rqa_result[[i]])
  }
  
  # Test 3
  ts = readRDS("../testdata/rqa_ts_3.RDS")  
  handmade_result = readRDS("../testdata/rqa_results_3.RDS")
  rqa_result = rqa(
    time.series = ts,
    embedding.dim = 2,
    time.lag = 1,
    radius = 1.2,
    lmin = 2,
    do.plot = FALSE,
    distanceToBorder = 1
  )
  
  for (i in names(handmade_result)) {
    expect_equal(handmade_result[[i]], rqa_result[[i]])
  }
})
