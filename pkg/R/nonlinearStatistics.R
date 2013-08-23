# #### 
# #### Computes the Takens' estimator.
# ####   
# #### details
# #### The Takens' estimator is an alternative to determining a finite correlation dimension suggested by Takens.
# #### Takens' estimator provides a maximum likelihood estimate and it is also used as statistic for nonlinearity
# #### testing.
# #### references Kantz, H.  and Schreiber, T.: Nonlinear Time series Analysis (Cambridge university press)
# #### Theiler, J. Lacunarity in a best estimator of fractal dimension. Phys. Lett. A,135, 195 (1988)
# ####
# #### author Constantino A. Garcia
# 
# takens.estimator = function(time.series, min.embedding.dim = 2, max.embedding.dim = 5, time.lag = 1,
#                             min.radius, max.radius, n.points.radius, do.plot = TRUE, 
#                             theiler.window = 0, number.boxes = NULL){
#   
#   
#     corr.dim = corrDim(time.series=time.series,min.embedding.dim= min.embedding.dim,
#                        max.embedding.dim=max.embedding.dim, time.lag = time.lag,
#                         min.radius = min.radius, max.radius = max.radius,
#                         n.points.radius = n.points.radius, do.plot=do.plot,
#                         theiler.window = theiler.dist, number.boxes =number.boxes)
# 
#     # vector that will store the takens' estimator for each dimension
#     takens.estimator = vector(mode = "numeric", length = nrow(corr.dim) )
#     names(takens.estimator) = min.embedding.dim:max.embedding.dim
#     # the radius axis. We have to revert it since it's sorted in decreasing order
#     radius.vector =  rev(as.numeric(dimnames(corr.dim)[[2]]))
#     # for each dimension of the correlation dimension matrix, we will calculate
#     # the takens' estimator
#     position = 1
#     for (dimension in min.embedding.dim:max.embedding.dim){
#       current.dimension = as.character(dimension)
#       integrand = rev(corr.dim[current.dimension,])/radius.vector
#       integral = trapezoidalRule(c(0,radius.vector), c(0,integrand))
#       lc=log(corr.dim[current.dimension,])
#       le=log(rev(radius.vector))
#       cat("*\n")
#       takens.estimator[[current.dimension]] = corr.dim[current.dimension,1]/integral
#     }
#     print(takens.estimator)
#     return (list(takens.estimator=takens.estimator,correlation.matrix = corr.dim))
# }




timeReversibility <- function(time.series,tau){
  len.ts = length(time.series)
  statistic = c()
  for (i in seq_along(tau)){
    statistic[[i]] = mean( (time.series[(tau[[i]]+1):len.ts] - time.series[1:(len.ts-tau[[i]])])^3)  
  }
  return(statistic)
}