test.poincare.section= function(){
  Ts=0.001
  r=rossler(a = 0.2, b = 0.2, c = 5.7, start=c(-2, -10, 0.2) ,time=seq(0,300,by = Ts), do.plot=FALSE)
  time.series=cbind(r$x,r$y,r$z)
  # calculate poincare sections
  cat("\n",getwd(),"\n")
  pm=poincareMap(takens = time.series,  normal.hiperplane.vector = c(0,1,0),  hiperplane.point=c(0,0,0) )
  pm2=poincareMap(takens = time.series, normal.hiperplane.vector = c(-1,1,0),  hiperplane.point=c(0,0,0) )
  pm3=poincareMap(takens = time.series, normal.hiperplane.vector = c(-1,0,0),  hiperplane.point=c(0,0,0) )
  pm4=poincareMap(takens = time.series, normal.hiperplane.vector = c(1,1,0),  hiperplane.point=c(0,0,0) )
  # check against stored results
  load(file="auxiliar_test_files/poincare.test.1")
  load(file="auxiliar_test_files/poincare.test.2")
  load(file="auxiliar_test_files/poincare.test.3")
  load(file="auxiliar_test_files/poincare.test.4")
  
  checkEquals(poincare.test.1,pm)
  checkEquals(poincare.test.2,pm2)
  checkEquals(poincare.test.3,pm3)
  checkEquals(poincare.test.4,pm4)
  
}
  