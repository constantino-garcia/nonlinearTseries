test.correlation.dimension = function(){
  ###################################################### henon   ###################################################################
  set.seed(1)
  h=henon(n.sample=  3000,n.transient= 100, a = 1.4, b = 0.3, 
          start = c(0.73954883, 0.04772637), do.plot = FALSE)
  
  ts=h$x
  mmin=2
  mmax=5
  time.lag=1
  rmin=10^-2.35
  rmax=10^-2
  np=20
  theiler.window=5
  x=corrDim(time.series=ts,min.embedding.dim=mmin,max.embedding.dim=mmax,time.lag=time.lag,min.radius=rmin,max.radius=rmax,n.points.radius=np,do.plot=FALSE,theiler.window=theiler.window,number.boxes=100)
  # cat("expected: 1.26  estimate: ",estimate(x),"\n")
  checkEqualsNumeric(1.250462,estimate(x),tolerance=10^-2)
  
  rmin=0.001
  rmax=0.01
  h=henon(n.sample=  3000,n.transient= 100, a = 1.2, b = 0.3, 
          start = c(1,1), do.plot = FALSE)
  
  ts=h$x
  x=corrDim(time.series=ts,min.embedding.dim=mmin,max.embedding.dim=mmax,time.lag=time.lag,min.radius=rmin,max.radius=rmax,n.points.radius=np,do.plot=FALSE,theiler.window=theiler.window,number.boxes=100)
  # cat("expected: 1.200  estimate: ",estimate(x),"\n")
  checkEqualsNumeric(1.211371,estimate(x),tolerance=10^-2)
  
  ############################################################# lorenz ################################################################
  
  l=lorenz(sigma=10, rho = 28, beta =8/3, start = c(-10, -11, 47), time =  seq(0, 70, by = 0.01), do.plot = FALSE)
  ts=l$x
  mmin=3
  mmax=6
  time.lag=10
  rmin=10^-0.6
  rmax=1
  np=5
  theiler.window=100
  x=corrDim(time.series=ts,min.embedding.dim=mmin,max.embedding.dim=mmax,time.lag=time.lag,min.radius=rmin,max.radius=rmax,n.points.radius=np,do.plot=FALSE,theiler.window=theiler.window,number.boxes=100)
  #cat("expected: 2.05  estimate: ",estimate(x),"\n")
  checkEqualsNumeric(2.095474,estimate(x),tolerance=6*10^-2)
  ######################################################## rossler ###################################################################
  r=rossler(a=0.15,b=0.2,c=10,start=c(0,0,0), time=seq(0,900,0.1),do.plot=FALSE)
  ts=r$x
  mmin=2
  mmax=5
  time.lag=12
  rmin=10^-0.55
  rmax=10^-0.2
  np=10
  theiler.window=100
  x=corrDim(time.series=ts,min.embedding.dim=mmin,max.embedding.dim=mmax,time.lag=time.lag,min.radius=rmin,max.radius=rmax,n.points.radius=np,do.plot=FALSE,theiler.window=theiler.window,number.boxes=100)
  #cat("expected: 1.26  estimate: ",estimate(x),"\n")
  checkEqualsNumeric(1.280408,estimate(x),tolerance=3*10^-2)
    
  ########################################################## logistic map #############################################################
  logmap=logisticMap(r=3.5699456,n.sample=5000,n.transient=500,do.plot=FALSE)
  ts=logmap
  x=corrDim(time.series=ts,min.embedding.dim=2,max.embedding.dim=4,time.lag=1,min.radius=10^-5,max.radius=10^-3.5,n.points.radius=10,do.plot=FALSE,theiler.window=100,number.boxes=100)
  #cat("expected: 0.538  estimate: ",estimate(x),"\n")
  checkEqualsNumeric(0.5057,estimate(x),tolerance=3*10^-2)
  
}