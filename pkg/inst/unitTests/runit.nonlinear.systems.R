test.nonlinear.systems=function(){
  
  # load files checked by visual inspection
  load(file="auxiliar_test_files/henon")
  load(file="auxiliar_test_files/logistic.map")
  load(file="auxiliar_test_files/lorenz")
  load(file="auxiliar_test_files/rossler")
  load(file="auxiliar_test_files/ikeda")
  load(file="auxiliar_test_files/clifford")
  load(file="auxiliar_test_files/sinai")
  load(file="auxiliar_test_files/gauss")
  
  cat("\n")
  cat("Generating Henon map\n")
  h.test=henon(n.sample = 1000, n.transient=10,do.plot=FALSE,start=c(-0.006423277,-0.473545134))
  checkEquals(h,h.test)
  
  cat("Generating logistic map\n")
  lm.test=logisticMap( n.sample=1000,n.transient=10,do.plot=FALSE,start=0.0234)
  checkEquals(lm,lm.test)
  
  cat("Generating Lorenz attractor\n")
  lor.test=lorenz(time=seq(0,30,by = 0.01),do.plot=FALSE)
  checkEquals(lor,lor.test)
  
  cat("Generating Rossler attractor\n")
  r.test=rossler(time=seq(0,30,by = 0.01),do.plot=FALSE,)
  checkEquals(r,r.test)
  
  cat("Generating Ikeda map\n")
  ik.test=ikedaMap(n.sample = 2000, n.transient=100,do.plot=FALSE,start=c(0.2345,-0.2438))
  checkEquals(ik,ik.test)
  
  cat("Generating Clifford map\n")
  cm.test=cliffordMap(n.sample = 2000, n.transient=100,do.plot=FALSE,start=c(0.3046741,0.4720825))
  checkEquals(cm,cm.test)
  
  cat("Generating Sinai map\n")
  sm.test=sinaiMap(n.sample = 5000, n.transient=100,do.plot=FALSE,start=c(0.2568524,0.3846888),a=0.3)
  checkEquals(sm,sm.test)
  
  cat("Generating Gauss map\n")
  gm.test=gaussMap(n.sample = 2000, n.transient=100,do.plot=FALSE,b=-0.5,start=0.3)
  checkEquals(gm,gm.test)
  
}