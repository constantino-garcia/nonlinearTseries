test.rqa = function(){
  load("auxiliar_test_files/completeTest1");load("auxiliar_test_files/paramsCT1")
  aa=rqa(time.series = ts, embedding.dim=2, time.lag=1,radius=1.2,lmin=2,do.plot=FALSE,distanceToBorder=2)
  n=names(a)
  for (i in names(a)){
    checkEquals(a[[i]],aa[[i]])
  }
  ###################################### test 3 ##################################
  load("auxiliar_test_files/paramsCT2")
  aa=rqa(time.series = ts, embedding.dim=2, time.lag=1,radius=1.2,lmin=3,vmin=4,do.plot=FALSE,distanceToBorder=3)
  
  for (i in names(a)){
    checkEquals(a[[i]],aa[[i]])
  }
  
  
  ###################################### test 3 ##################################
  load("auxiliar_test_files/completeTest3");load("auxiliar_test_files/paramsCT3")
  aa=rqa(time.series = ts, embedding.dim=2, time.lag=1,radius=1.2,lmin=2,do.plot=FALSE,distanceToBorder=1)
  
  for (i in names(a)){
    checkEquals(a[[i]],aa[[i]])
  }
  
 
  ts=rnorm(200);tak  = buildTakens(ts,2,1)
  mat = as(as.matrix(dist(tak,method="max"))<0.1,"sparseMatrix");image(mat)
  recurrencePlot(NULL,ts,embedding.dim=2,1,0.1)
  
}