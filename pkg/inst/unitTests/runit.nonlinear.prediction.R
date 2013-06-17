test.nonlinear.prediction = function(){
  h=henon(n.sample=5000,start=c(0.324,-0.8233))

  ################################# basic test ##################################
  m=2
  lag=1
  start=2600
  predictionStep=1
  epsi=0.1
  ts=h$x[1:2500]
  predic=nonLinearPrediction(time.series=ts,embedding.dim=2,
                             time.lag=lag,
                             prediction.step=predictionStep,radius=epsi,radius.increment=epsi/2)
  
  cat("real value: ",h$x[2500+predictionStep],"n")
  cat("prediction: ",predic,"\n")
  
}