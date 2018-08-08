library(nonlinearTseries)
context("ACF-based time-lag estimation")

test_that("the time lag for the Takens' vectors are properly estimated", {
  skip_on_cran()
  set.seed(1)
  Ts = 0.01
  time = seq(0, 100, Ts)
  time.series = cos(2 * pi * 1 * time)
  tl = timeLag(
    time.series,
    technique = "acf",
    selection.method = "first.zero",
    value = 0,
    lag.max = NULL,
    do.plot = FALSE
  )
  expect_equal(25, tl)
  
  time.series = exp(-time)
  tl = timeLag(
    time.series,
    technique = "acf",
    selection.method = "first.e.decay",
    value = 0,
    lag.max = 120,
    do.plot = FALSE
  )
  expect_equal(100, tl)
  time.series = arima.sim(n = 5000, model = list(ar = 0.9, -0.9))
  tl = timeLag(
    time.series,
    technique = "ami",
    n.partitions = NULL,
    selection.method = "first.e.decay",
    value = 0,
    lag.max = 120,
    do.plot = FALSE
  )
  expect_equal(1, tl)
  
  
  time.series = time
  tl = timeLag(
    time.series,
    technique = "acf",
    selection.method = "first.minimum",
    value = 0,
    lag.max = 8000,
    do.plot = FALSE
  )
  expect_equal(7072, tl)
  tl = timeLag(
    time.series,
    technique = "ami",
    n.partitions = NULL,
    selection.method = "first.minimum",
    value = 0,
    lag.max = 400,
    do.plot = FALSE
  )
  expect_equal(246, tl)
  
  
  
  tl = timeLag(
    time.series,
    technique = "acf",
    selection.method = "first.value",
    value = 0.8,
    lag.max = 8000,
    do.plot = FALSE
  )
  expect_equal(669, tl)
  tl = timeLag(
    time.series,
    technique = "ami",
    n.partitions = NULL,
    selection.method = "first.value",
    value = 2.5,
    lag.max = 400,
    do.plot = FALSE
  )
  expect_equal(120, tl)

})

test_that("exceptions are properly generated", {
  time.series = seq(0, 100, 0.01)
  # function does not cross 0.8
  expect_error(
    timeLag(
      time.series,
      technique = "acf",
      selection.method = "first.value",
      value = 0.8,
      lag.max = 10,
      do.plot = FALSE
    )
  )
  # function does not have a minimum
  expect_error(
    timeLag(
      time.series,
      technique = "acf",
      selection.method = "first.minimum",
      lag.max = 10,
      do.plot = FALSE
    )
  )
})
