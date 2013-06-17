test.information.dimension = function(){
  q=2
  ###################################################### henon   ###################################################################
  h=henon(n.sample=  3000,n.transient= 100, a = 1.4, b = 0.3, 
          start =  c(0.8253681, 0.6955566), do.plot = FALSE)
  
  ts=h$x

  leps = infDim(ts,2,1, min.fixed.mass=0.04, max.fixed.mass=0.2,
                              number.fixed.mass.points=100, radius =0.001, 
                              increasing.radius.factor = sqrt(2), number.boxes=100, 
                              number.reference.vectors=100, theiler.window = 10, 
                              kMax = 100,do.plot=FALSE)
  #cat("Henon---> expected: 1.24    predicted: ",estimate(leps),"\n")
  checkEqualsNumeric(1.178795,estimate(leps,do.plot=FALSE),tolerance=10^-2)
  
  
  
  ##################################### rossler #################################
  r=rossler(a = 0.2, b = 0.2, c = 5.7, start=c(-2, -10, 0.2)
            ,time=seq(0,60,by = 0.01), do.plot=FALSE)
  ts=r$x
  leps = infDim(ts,5,100, min.fixed.mass=0.001, max.fixed.mass=0.1,
                              number.fixed.mass.points=100, radius =0.1, 
                              increasing.radius.factor = sqrt(2), number.boxes=100, 
                              number.reference.vectors=100, theiler.window = 200, 
                              kMax = 100,do.plot=FALSE)
  #cat("Rossler--->expected: 2.92    predicted: ",estimate(leps),"\n")
  checkEqualsNumeric(3.00,estimate(leps,do.plot=FALSE),tolerance=5*10^-2)
  ################################# sinai map #########################################
  s=sinaiMap(a=0.3,n.sample=5000,start=c(0.23489,0.8923),do.plot=FALSE)
  ts=s$x

  leps = infDim(ts,2,1, min.fixed.mass=0.01, max.fixed.mass=0.03,
                              number.fixed.mass.points=1000, radius =0.1, 
                              increasing.radius.factor = sqrt(2), number.boxes=100, 
                              number.reference.vectors=100, theiler.window = 200, 
                              kMax = 100,do.plot=FALSE)
  checkEqualsNumeric(1.68,estimate(leps,do.plot=FALSE),tolerance=5*10^-2)
  #cat("Rossler--->expected: 1.72-1.77    predicted: ",estimate(leps),"\n")
  
}