

vectorizePar = function(par, N, default=1:N){
  if (is.null(par)) par=default
  if (length(par) < N) par=rep(par,length.out=N)
  par
}
                        
