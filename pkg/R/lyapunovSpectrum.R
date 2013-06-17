# falta cambiar los parametros por defecto metidos a pelo por otros mas adecuiados!!!!"!!
################################################### Biblio ########################################################
# based on Eckmann, Ruelle: ergodic  theory of chaos and sano-sawada-Measurement of the Lyapunov Spectrum from a Chaotic Time Series
###########################################################################################################
lyapunovSpectrum=function(takens,radius,radius.increment,number.boxes=NULL,nneighs=5,n.trajectories=50,
                          max.distance.from.reference=10,sampling.period=1, do.plot=TRUE){
  #estimates number.boxes if it has not been specified
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  # auxiliar variables
  ntakens = nrow(takens)
  # one file for each lyapunov exponent. The k-th colum represents the value
  # of every lyapunov exponent after k steps from the reference point
  lyapunov = matrix(0,nrow=ncol(takens),ncol=max.distance.from.reference)
  # last takens' vector that we can use as a initial point for a trajectory
  last.takens = (ntakens-max.distance.from.reference)
  trajectories.used =  0
  # check
  if (n.trajectories > last.takens) stop("number of trajectories greater than available takens' vectors\n")
  # compute the start points for the n.trajectories
  delta = (last.takens+1)/n.trajectories
  start.points = seq(1,last.takens,by=delta)
  # find at least nneighs neighbours for each takens' vector
  unnnn = proc.time()
  cat("searching neighs\n")
  neighbours = find.at.least.N.neighboursBIS(takens,radius,radius.increment,nneighs,number.boxes)
  print(proc.time()-unnnn)
  cat("jacobians\n")
  unnnn = proc.time()
  # estimate the linearized flow maps for every point in a trajectory (except the last point in each
  # trajectory, that we won't use to predict because there is no next step)
  auxiliar.next.points = 0:(max.distance.from.reference-1)
  trajectory.points = start.points[[1]] + auxiliar.next.points
  n.start.points = length(start.points)
  for (j in 2:n.start.points){
    trajectory.points = union(trajectory.points,start.points[[j]] + auxiliar.next.points)
  }
  linearized.flow.maps = estimateLinearizedFlowMap(takens,trajectory.points,neighbours)
  print(proc.time()-unnnn)
  cat("averaging\n")
  unnnn = proc.time()
  for (sp in start.points){
    # compute the lyapunov spectrum using the local jacobians in the trajectory
    # that starts in the actual starting point
    aux.lyapunov = lyapunovFromTrajectory(takens,sp,linearized.flow.maps,max.distance.from.reference,sampling.period)
    # avoid infinites due to 0s in some Ri matrix
    if (!any(is.infinite(aux.lyapunov))) {
      # store to average
      lyapunov = lyapunov + aux.lyapunov
      trajectories.used = trajectories.used + 1
    }
    
  }
  
  # compute the average lyapunov exponents
  if (trajectories.used > 0) lyapunov = lyapunov / trajectories.used
  # plotting
  if (do.plot){
    n.col.lyap = nrow(lyapunov)
    axis = (1:max.distance.from.reference)*sampling.period
    plot(axis,lyapunov[1,],col=1,ylim=range(lyapunov),xlab="Time (t)", 
         ylab=expression("Lyapunov exponents ("*lambda*")"),
         main=expression("Time evolution of Lyapunov exponents ("*lambda*"(t))"),cex=0.4,type="b")
    for (i in 2:n.col.lyap){
      lines(axis,lyapunov[i,],col=i,cex=0.4,type="b")
    }
  }
  print(proc.time()-unnnn)
  # returning lyapunov exponents
  return (lyapunov)
}

estimateLinearizedFlowMap=function(takens,trajectory.points,neighbours){
  embedding.dim=ncol(takens)
  n.takens = nrow(takens)
  # arrays to store the linearized flow maps
  linearized.flow.maps=array(0,dim=c(n.takens,embedding.dim,embedding.dim))
  
  for (point in trajectory.points){
    # prepare vectors of close neighbours in present and in the next step
    aux.neighbours = neighbours[[point]]
    neighs.found = length(aux.neighbours)    
    if (neighs.found > 0){
      # avoid get an 'index out of bounds' error (we will search over the next position, thus
      # we will avoid the last position of the takens' vectors array)
      valid.neighbours = aux.neighbours[aux.neighbours < nrow(takens)]
      n.valid.neighbours=length(valid.neighbours)
      # compute de delta.vectors computations
      delta.vector = matrix(0,nrow=n.valid.neighbours,ncol=embedding.dim )   
      evolution.delta.vector = matrix(0,nrow=n.valid.neighbours,ncol=embedding.dim )   
      for (it.delta in 1:n.valid.neighbours){
        delta.vector[it.delta,]=takens[valid.neighbours[[it.delta]],]-takens[point,]
        evolution.delta.vector[it.delta,]=takens[valid.neighbours[[it.delta]]+1,]-takens[point+1,]
      }
    }else{ 
      # we have found and isolated takens' vector (with no neighbours!!!)
      stop("There exist a takens' vector with no neighbours!!! Quitting!!\n")
    }
    
    ######################### estimate the local linearized flow map ##################
    ######################### see sano and sawada ###########################
    # compute the autocovariance matrix of the delta.vectors to estimate the  linearized  local flow map 
    autocov = t(delta.vector)%*%delta.vector
    # compute the cross-covariance matrix of the delta.vectors and the evolution.delta.vector
    # to estimate the  linearized local flow map
    cross.cov = t(evolution.delta.vector)%*%delta.vector
    # use the pseudoinverse to avoid problems with singular autocov
    linearized.flow.maps[point,,] = cross.cov %*% pseudoinverse(autocov)
    # increase the number of the reference.vector
  }
  return(linearized.flow.maps)
  
}


lyapunovFromTrajectory = function(takens,number.of.reference.vector,linearized.flow.maps,
                                  max.distance.from.reference,sampling.period){
  embedding.dim=ncol(takens)
  # arrays to store Q and R matrix from the QR factorization
  Q=array(0,dim=c(max.distance.from.reference+1,embedding.dim,embedding.dim))
  R=array(0,dim=c(max.distance.from.reference,embedding.dim,embedding.dim))
  lyapunov.exponents=matrix(0,nrow=embedding.dim,ncol=max.distance.from.reference)
  Q[1,,]=diag(embedding.dim)
  
  step = 1
  trajectory = number.of.reference.vector:(number.of.reference.vector + (max.distance.from.reference-1))
  for (point in trajectory){
    
    local.linearized.flow.map=linearized.flow.maps[point,,]
    ######################### estimate the lyapunov spectrum ##################
    ######################### see Eckmann and Ruelle ###########################
    # compute the evolution of an orthonomal set of vectors (Q matrix has an 
    # orthonormal vector per colum)
    evolution = local.linearized.flow.map %*% Q[step,,]
    # perform a qr factorization to obtain a new set of orthonormal vectors
    qr.factorization = qr(evolution)
    Q[step+1,,] = qr.Q(qr.factorization)
    R[step,,] = qr.R(qr.factorization)
    # transform the R matrix so it has positive diagonal components (in that way
    # that the factorization is unique and also, to avoid problems when applying
    # the logarithm to the diagonal components)
    S = diag(sign(diag(R[step,,])))
    R[step,,] = S %*% R[step,,]
    Q[step+1,,] = Q[step+1,,] %*% S
    
    
    if (step == 1) {
      lyapunov.exponents[,step] = log( diag(R[step,,]) )
    }else{
      lyapunov.exponents[,step] = lyapunov.exponents[,step-1] + log( diag(R[step,,]) )
    }
    #increase step, the index for storing Q and R matrix
    step = step +1
  }
  
  # average each column
  for (col in 2:max.distance.from.reference){
    lyapunov.exponents[,col] = lyapunov.exponents[,col]/(col*sampling.period)
  }
  
  return(lyapunov.exponents)
}
