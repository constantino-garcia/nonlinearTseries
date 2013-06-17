test.if.runs = function(){
  
  # there is no test yet for these functions...
  # at least we will check that they run properly
  
  
  N=3000
  sd=0.25#noise power
  
  h=henon(n.sample=N, start=c(-0.006423277,-0.473545134))
  noise=rnorm(N,sd=sd)
  h$x=h$x+noise
  h$y=h$y+noise
  m=2
  lag=1
  predictionStep=1
  epsi=0.1
  denx=nonLinearNoiseReduction(h$x,1*sd,embedding.dim=m)
  deny=nonLinearNoiseReduction(h$y,1*sd,embedding.dim=m)
  
  
  
}