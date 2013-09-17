test.timeLag = function(){
  Ts = 0.01
  time = seq(0,100,Ts)  
  time.series = cos(2*pi*1*time)
  tl = timeLag(time.series, method="first.zero",value = 0,lag.max=NULL,do.plot=TRUE)  
  checkEquals(25,tl)
  
  time.series = exp(-time)
  tl = timeLag(time.series, method="first.e.decay",value = 0,lag.max=NULL,do.plot=TRUE)  
  checkEquals(100,tl)
  
  time.series = time
  tl = timeLag(time.series, method="first.minimum",value = 0,lag.max=NULL,do.plot=TRUE)  
  checkEquals(7072,tl)
  
  
  tl = timeLag(time.series, method="first.value",value = 0.8,lag.max=NULL,do.plot=TRUE)  
  checkEquals(669,tl)
  
  cat("\n...checking if exceptions are given...\n")
  checkException(timeLag(time.series, method="first.value",value = 0.8,lag.max=10,do.plot=TRUE)  )
  checkException(timeLag(time.series, method="first.minimum",lag.max=10,do.plot=TRUE)  )
}