binScore2 <- function(cM, rM){
  comp <- (cM[1,1]>rM[1,1])&&(cM[1,2]<rM[1,2])&&(cM[1,3]<=rM[1,3]) 
  (cM[2,1]>rM[2,1])&&(cM[2,2]<rM[2,2])&&(cM[2,3]< rM[2,3])
  (cM[3,1]>rM[3,1])&&(cM[3,2]<rM[3,2])&&(cM[3,3]<=rM[3,3])
  return(comp)
}

toReqMat <- function (aCountsMat, aReqPercentMat){
  return(round(aReqPercentMat*sum(aCountsMat)/100,0))
}

reqPercentages <- matrix (c(15, 5, 0, 10, 5, 5, 10, 10, 15), nrow=3, byrow=TRUE)

(countsRnd    <- matrix(floor(runif(9)*10)+1, nrow=3, ncol=3))
(reqRnd <- toReqMat (countsRnd, reqPercentages)) 
binScore2(countsRnd, reqRnd)

(countsTrueL1 <- matrix (c(20, 3, 0, 10, 2, 2, 20, 10, 20), nrow=3, byrow=TRUE))
(reqTrueL <- toReqMat (countsTrueL1, reqPercentages))
binScore2(countsTrueL1, reqTrueL)
