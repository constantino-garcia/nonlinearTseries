test.dfa = function(){
  set.seed(1)
  noise=rnorm(6000)
  brown_noise=function(len=100){
    br=c()
    x=rnorm(len)
    for (i in 1:len){
      br[[i]]=sum(x[1:i])
    }
    return(br)
  }
  
  bnoise=brown_noise(6000)
  #load stored pink noise
  load(file="auxiliar_test_files/pink.noise.ts")
  ####### DFA
  sm=10;sM=1000
  miDFA=dfa(time.series=noise,npoints=10,window.size.range=c(sm,sM),do.plot=FALSE);
  cat("teórico: 0.5---Estimado: ", estimate(miDFA),"\n")
  checkTrue(abs(0.50-estimate(miDFA))/0.50 < 0.1)
  #cat("teórico: 0.5---Estimado: ",miDFA$alpha1,"\n")
  
  miDFA=dfa(time.series=pnoise,npoints=10,window.size.range=c(sm,sM),do.plot=FALSE);
  #cat("teórico: 1---Estimado: ",miDFA$alpha1,"\n")
  checkTrue(abs(1-estimate(miDFA))/1< 0.1)
 
  
  miDFA=dfa(time.series=bnoise,npoints=10,window.size.range=c(sm,sM),do.plot=FALSE);
  checkTrue(abs(1.50-estimate(miDFA))/1.50 < 0.1)
  #cat("teórico: 1.5---Estimado: ",miDFA$alpha1,"\n")
  

}