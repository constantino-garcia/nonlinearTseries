library(nonlinearTseries)
context("RQA")

test_that("test rqa", {
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
