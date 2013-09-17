test.neighbour.search=function(){
  eps=0.4
  cat("\n")
  for (k in c(2,4,5,10)){
    load(file=paste("auxiliar_test_files/takensVector_",k,sep=""))
    name = paste("auxiliar_test_files/neighbourSearchSolution_",k,sep="")
    load(file=name)
    cat("searching for neighbours in a ",k,"-dimensional embedding\n")
    sc=findAllNeighbours(data,eps)
    for (j in 1:nrow(data)){
      checkEquals(length(s[[j]]),length(sc[[j]]))
      if (length(s[[j]])!=0){
        checkEquals(s[[j]], sort(sc[[j]]) )
      }
    }
  }
  
}
