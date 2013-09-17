test.takens.estimator = function(){
  
  ########################### Henon ##########################################
#   cat("\nTakens estimator for the Henon attractor...\n")
#   h=henon(n.sample=  3000,n.transient= 100, a = 1.4, b = 0.3, 
#           start = c(0.73954883, 0.04772637), do.plot = FALSE)
#   
#   ts=h$x
#   mmin=2
#   mmax=5
#   time.lag=1
#   rmin=10^-3
#   rmax=10^-2
#   np=1000
#   tdist=5
#   x=takens.estimator(ts, min.embedding.dim = mmin, max.embedding.dim = mmax, time.lag = time.lag,
#                      rmin, rmax, np, do.plot = TRUE, 
#                      theiler.window = tdist, number.boxes = NULL)
#   cat("expected: 1.26  estimate: ",median(x$takens.estimator),"\n")
#   x=takens.estimator(ts, min.embedding.dim = mmin, max.embedding.dim = mmax/2, time.lag = time.lag,
#                      rmin, rmax, np, do.plot = TRUE, 
#                      theiler.window = tdist, number.boxes = NULL)
#   cat("expected: 1.26  estimate: ",median(x$takens.estimator),"\n")
  #checkEqualsNumeric(1.250462,x,tolerance=10^-2)
  
#   cat("Takens estimator for the Henon attractor...\n")
#   min=0.001
#   rmax=0.01
#   h=henon(n.sample=  3000,n.transient= 100, a = 1.2, b = 0.3, 
#           start = c(1,1), do.plot = FALSE)
#   
#   ts=h$x
#   x=takens.estimator(ts, min.embedding.dim = mmin, max.embedding.dim = mmax, time.lag = time.lag,
#                      rmin, rmax, np, do.plot = TRUE, 
#                      theiler.window = tdist, number.boxes = NULL)
#   cat("expected: 1.20  estimate: ",median(x$takens.estimator),"\n")
#   #checkEqualsNumeric(1.211371,estimate(x),tolerance=10^-2)
#   
#   ########################################################### lorenz ################################################################
#   cat("Takens estimator for the Lorenz attractor...\n")
#   l=lorenz(sigma=10, rho = 28, beta =8/3, start = c(-10, -11, 47), time =  seq(0, 70, by = 0.01), do.plot = FALSE)
#   ts=l$x
#   mmin=3
#   mmax=6
#   time.lag=10
#   rmin=10^-0.6
#   rmax=1
#   np=5
#   tdist=100
#   x=takens.estimator(ts, min.embedding.dim = mmin, max.embedding.dim = mmax, time.lag = time.lag,
#                      rmin, rmax, np, do.plot = TRUE, 
#                      theiler.window = tdist, number.boxes = NULL)
#   cat("expected: 2.05  estimate: ",median(x$takens.estimator),"\n")
#   #checkEqualsNumeric(2.095474,estimate(x),tolerance=6*10^-2)
#   ######################################################## rossler ###################################################################
#   cat("Rossler estimator for the Lorenz attractor...\n")
#   r=rossler(a=0.15,b=0.2,c=10,start=c(0,0,0), time=seq(0,900,0.1),do.plot=FALSE)
#   ts=r$x
#   mmin=2
#   mmax=5
#   time.lag=12
#   rmin=10^-0.55
#   rmax=10^-0.2
#   np=10
#   tdist=100
#   x=takens.estimator(ts, min.embedding.dim = mmin, max.embedding.dim = mmax, time.lag = time.lag,
#                      rmin, rmax, np, do.plot = TRUE, 
#                      theiler.window = tdist, number.boxes = NULL)
#   cat("expected: 1.26  estimate: ",median(x$takens.estimator),"\n")
#   print(x)
#   #checkEqualsNumeric(1.280408,estimate(x),tolerance=3*10^-2)
#   
#   ######################################################### logistic map #############################################################
#   mmin=2
#   mmax=5
#   time.lag=12
#   rmin=10^-0.55
#   rmax=10^-0.2
#   np=10
#   tdist=100
#   cat("Logistic map estimator for the Lorenz attractor...\n")
#   logmap=logisticMap(r=3.5699456,n.sample=5500,n.transient=500,do.plot=FALSE)
#   ts=logmap
#   x=takens.estimator(ts, min.embedding.dim = mmin, max.embedding.dim = mmax, time.lag = time.lag,
#                      rmin, rmax, np, do.plot = TRUE, 
#                      theiler.window = tdist, number.boxes = NULL)
#   cat("expected: 0.538  estimate: ",median(x$takens.estimator),"\n")
#   #checkEqualsNumeric(0.5057,estimate(x),tolerance=3*10^-2)
#   
}
