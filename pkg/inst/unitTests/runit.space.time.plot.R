test.space.time.plot = function(){
  
  tak=buildTakens(cos(2*pi*0.01*0:5000),2,1)
  stp.test=spaceTimePlot( takens=tak,number.time.steps=400,do.plot=FALSE)
  load(file="auxiliar_test_files/cosine_stp")
  checkEquals(stp[[1]],stp.test$stp.matrix[[1]])
  checkEquals(stp[[2]],stp.test$stp.matrix[[2]])
  checkEquals(stp[[3]],stp.test$stp.matrix[[3]])
  
  tak=buildTakens(sin(2*pi*0.01*(0:5000)),2,1)
  stp.test=spaceTimePlot(takens=tak,number.time.steps=400,do.plot=FALSE)
  load(file="auxiliar_test_files/sine_stp")
  checkEquals(stp[[1]],stp.test$stp.matrix[[1]])
  checkEquals(stp[[2]],stp.test$stp.matrix[[2]])
  checkEquals(stp[[3]],stp.test$stp.matrix[[3]])
}