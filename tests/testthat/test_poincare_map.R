library(nonlinearTseries)
context("Poincare maps")

test_that("test Poincare maps", {
  skip_on_cran()
  set.seed(1)
  Ts=0.001
  r = rossler(a = 0.2, b = 0.2, w = 5.7,
    start = c(-2, -10, 0.2), time = seq(0, 300, by = Ts),
    do.plot = FALSE
  )
  time.series = cbind(r$x, r$y, r$z)
  # calculate poincare sections
  pm = poincareMap(
    takens = time.series,
    normal.hiperplane.vector = c(0, 1, 0),
    hiperplane.point = c(0, 0, 0)
  )
  pm2 = poincareMap(
    takens = time.series,
    normal.hiperplane.vector = c(-1, 1, 0),
    hiperplane.point = c(0, 0, 0)
  )
  pm3 = poincareMap(
    takens = time.series,
    normal.hiperplane.vector = c(-1, 0, 0),
    hiperplane.point = c(0, 0, 0)
  )
  pm4 = poincareMap(
    takens = time.series,
    normal.hiperplane.vector = c(1, 1, 0),
    hiperplane.point = c(0, 0, 0)
  )
  
  # check against stored results that have been visually validated
  load(file = "../testdata/poincare.test.1")
  load(file = "../testdata/poincare.test.2")
  load(file = "../testdata/poincare.test.3")
  load(file = "../testdata/poincare.test.4")
  
  expect_equal(poincare.test.1,pm[1:6])
  expect_equal(poincare.test.2,pm2[1:6])
  expect_equal(poincare.test.3,pm3[1:6])
  expect_equal(poincare.test.4,pm4[1:6])
  
})
