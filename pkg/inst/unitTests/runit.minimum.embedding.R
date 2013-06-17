test.minimum.embedding = function(){
  

h=henon(do.plot=FALSE,start=c(0.4,0.234),n.sample=3000)

l=lorenz(sigma=10, rho = 28, beta =8/3, start = c(-10, -11, 47), time =  seq(0, 30, by = 0.01), do.plot = FALSE)

ik=ikedaMap(n.sample=3000,n.transient=100,start=c(0.5341,0.278),do.plot=FALSE)

#brown_noise=rnorm(3000)
#for (i in 2:3000) brown_noise[[i]]=0.95*brown_noise[[i-1]]+brown_noise[[i]]

#results expected: 4,2,3,none
x=getEmbeddingDim(ik$x,time.lag=1,max.embedding.dim=6,theiler.window=10,threshold=0.9,do.plot=TRUE)
checkEquals(4,x)
x=getEmbeddingDim(h$x,time.lag=1,max.embedding.dim=6,theiler.window=10,threshold=0.9,do.plot=TRUE)
checkEquals(2,x)
x=getEmbeddingDim(l$x,time.lag=15,max.embedding.dim=6,theiler.window=50,threshold=0.9,do.plot=TRUE)
checkEquals(3,x)
#getEmbeddingDim(brown_noise,time.lag=1,max.embedding.dim=15,theiler.window=100,threshold=0.9)

}