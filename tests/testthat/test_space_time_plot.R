
library(nonlinearTseries)
context("space time plot port")

test_that("ported code yields same results", {
  nrepeat = 3
    for (i in seq_len(nrepeat)) {
      embeddingD = sample(2:10, 1)
      time.lag = sample(1:10, 1)
      
      number.boxes = sample(c(4, 10, 50, 100, 1000), 1)
      
      ts = arima.sim(model = list(ar = runif(1, -1, 1),
                                  ma = runif(1, -1 , 1)), 
                     1000)
      or = rcppSpaceTimePlot(time.series = ts, embedding.dim = embeddingD, time.lag = time.lag, do.plot=FALSE)
      ported = rcppSpaceTimePlot(time.series = ts, embedding.dim = embeddingD, time.lag = time.lag, do.plot=FALSE)
      
      expect_equal(or, ported)
      
    #   indx = which(matrix(is.na(or$stp.matrix),nrow = nrow(or$stp.matrix)),
    #                arr.ind = T)
    #   or$stp.matrix[indx]
    #   ported$stp.matrix[indx]
    } 
})
