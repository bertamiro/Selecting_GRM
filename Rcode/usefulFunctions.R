zeros <-function(x){which (x==0)}
discard <- function(x, where, howMany){
  return(length(zeros(x[where])) > howMany) 
}

discardA <- function(x, percentage){ 
  maxZeros <- ceiling (length(x)*percentage)
  return (discard(x,1:length(x),maxZeros))
}
# test
# exData <- matrix (0, nrow=10, ncol=10)
# for (i in 1:10){
#     for(j in 1:i)
#         exData[i,j]<-1
# }
# d1 <-apply(exData, 1, discardA, 0.66); length(d1); sum(d1);show(exData1 <- exData[!d1,])
# d1 <-apply(exData, 1, discardA, 0.5); length(d1); sum(d1);show(exData1 <- exData[!d1,])
# d1 <-apply(exData, 1, discardA, 0.33); length(d1); sum(d1);show(exData1 <- exData[!d1,])
# d1 <-apply(exData, 1, discardA, 0.2); length(d1); sum(d1);show(exData1 <- exData[!d1,])


countNAs <- function (X){
  sumaX <- apply(X, 1,sum); 
  # sum(is.na(sumaX))
  xAmbNAs <- X[is.na(sumaX),]
  NAs <- apply(xAmbNAs, 1, function(unVec){sum(is.na(unVec))})
}

DimMat <- function (y, title="", cols2show=5){
  show(dim(y))
  cat("\nFirst rows and columns\n")
  show(head(y[,1:cols2show]))
  cat("\nLast rows first columns\n")
  show(tail(y[,1:cols2show]))
#  cat("\nLast rows last columns\n")
#  show(tail(y[,(dim(y)[2]-cols2show+1):(dim(y)[2])]))
}