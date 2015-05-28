

vectorizePar = function(par, N, default=1:N){
  if (is.null(par)) par=default
  if (length(par) < N) par=rep(par,length.out=N)
  par
}
                        

propagateTakensAttr = function(x, takens){
  attr(x, "id") = attr(takens,"id")
  attr(x, "time.lag") = attr(takens, "time.lag")
  attr(x, "embedding.dim") = attr(takens, "embedding.dim")
  x
}
