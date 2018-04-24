require(Hmisc)
require(energy)
require(FactoMineR)

vecCorrs <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  return(unlist(list(rhoS=corrS$r[1,2],rhoP=corrP$r[1,2],
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

#' A row-wise correlation function calculator
#'
#' Given two matrices X (m,n), Y(m,n) this function computes (m) Pearson and Spearman correlation coefficients 
#' and their significance p-values for every pair of row vectors.
#' @param X First matrix
#' @param Y Second matrix. Must have the same dimensions as X
#' @keywords Correlation
#' @export
#' @examples 
#' (X <- round(matrix (rnorm(30)*10, ncol=6),1)) + 1:10
#' (Y <- round(X + matrix (rnorm(30)*10, ncol=6),1)) - 10:1
#' (rownames(X)=rownames(Y)=letters[1:nrow(X)])
#' (m1<-matCorrs(X,Y))
#' matCorrs
#' 
matCorrs <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  corrsList<- matrix(NA, nrow=nrow(X), ncol=4)
  colnames(corrsList) <- c("r (Sp)", "r (Pear)",  "p (Sp)","p (Pear)")
  for (i in 1:nrow(X)){
    corrs<- vecCorrs(X[i,],Y[i,])
    corrsList[i,] <- corrs
  }
  rownames(corrsList) <- rownames(X)
  return(corrsList)
}

distCorrs <- function (x, y){
  distCorr <- dcor(x, y, index=1.0)
  return(unlist(list(distCorr)))
}

matDistCorr <- function (X, Y){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  distCorrsList<- matrix(NA, nrow=nrow(X), ncol=1)
  colnames(distCorrsList) <- c("DistCor")
  for (i in 1:nrow(X)){
    DistCorr<- distCorrs(X[i,],Y[i,])
    distCorrsList[i,] <- DistCorr
  }
  rownames(distCorrsList) <- rownames(X)
  return(distCorrsList)
}

allCorrs  <- function (x, y){
  corrS <- rcorr(x, y, type="spearman")
  corrP <- rcorr(x, y, type="pearson")
  distCorr <- dcor(x, y, index=1.0)
  return(unlist(list(rhoS=corrS$r[1,2],  rhoP=corrP$r[1,2], distCorr,
                     pvalS=corrS$P[1,2], pvalP=corrP$P[1,2])))
}

sort1 <- function (X, col,DEC=TRUE, ...){
  return(X[sort.list(X[,col], decreasing=DEC), ])
}

matAllCorrs  <- function (X, Y, sortByCorrs = FALSE){
  if ((nrow(X)!=nrow(Y))||(ncol(X)!=ncol(Y))) stop('matrices dimensions do not match')
  corrsList<- matrix(NA, nrow=nrow(X), ncol=5)
  colnames(corrsList) <- c("r (Sp)", "r (Pear)", "distCor", "p (Sp)",  "p (Pear)")
  for (i in 1:nrow(X)){
    corrs<- allCorrs(X[i,],Y[i,])
    corrsList[i,] <- corrs
  }
  corrsList<- cbind(corrsList, adj.Spear.Pval= p.adjust(corrsList[,"p (Sp)"],"fdr"))
  corrsList<- cbind(corrsList, adj.Pear.Pval = p.adjust(corrsList[,"p (Pear)"],"fdr"))
  rownames(corrsList) <- rownames(X)
  if (sortByCorrs) corrsList <- sort1(corrsList,4, DEC=FALSE)
  return(corrsList)
}

naiveSelection <- function (X, Y, type="Spearman", 
                            adj=TRUE, pValCutoff=0.05, rCutoff=0,
                            sortByCorrs=FALSE){
  corsMat <- matAllCorrs (X, Y, sortByCorrs=sortByCorrs)
  if (type=="Spearman"){
      selected <- corsMat[,c("r (Sp)", "p (Sp)", "adj.Spear.Pval", "distCor")]
      if(adj){
        lShaped <-(selected[,"r (Sp)"]<rCutoff) & (selected[,"adj.Spear.Pval"] < pValCutoff)
      }else{
        lShaped <-(selected[,"r (Sp)"]<rCutoff) & (selected[,"p (Sp)"] < pValCutoff)
      }
  }else{
    selected <- corsMat[,c("r (Pear)", "p (Pear)", "adj.Pear.Pval", "distCor")]
    if(adj){
      lShaped <-(selected[,"r (Pear)"]< rCutoff) & (selected[,"adj.Pear.Pval"] < pValCutoff)
    }else{
      lShaped <-(selected[,"r (Pear)"]< rCutoff) & (selected[,"p (Pe)"] < pValCutoff)
    }
  }
  selected2 <- data.frame (selected, lShaped)
  colnames(selected2) <- c(colnames(selected), "SigNegCorr")
  return(selected2)
}

multivCorr <- function(X,Y){
  # stopifnot(require(energy))
  # stopifnot(require(FactoMineR))
  # RV1 <-coeffRV(X, Y)
  # cat("p-value : ", RV1$p.value,"\n")
  # dcor1<-dcor(X, Y)
  # cat("DistCorr: ", dcor1,"\n")
  coin1 <- cia(X, Y)
  cat("RV coeff: ", coin1$coinertia$RV,"\n")
}

Qnn<-function(x, q) # q: quantil en % (1-100)
{
  x<-as.numeric(x)
  Q<-quantile(x, q/100, na.rm=TRUE)
  return(Q)
}

#QInf i QSup determinen quin quartil volem calcular, metInf i metSup quin  el llindar dels quartils 
#que volem utilitzar (QInf<metInf, QSup>metSup)

Initialselection<-function(allData, metilData, corrData, 
                           QInf=25, metInf, QSup=75, metSup, 
                           Adjust=FALSE, pAdj)
{
  QSupDades<-apply(metilData,1,Qnn, QSup) ###de les metilacions
  QInfDades<-apply(metilData,1,Qnn, QInf)  ###de les metilacions
  pv<-corrData[,5]
  if (Adjust)   pv<-corrData[,7]
  selCond  <- corrData[ ,2] < 0 & 
    pv < pAdj & 
    QSupDades > metSup & 
    QInfDades < metInf
  Cors <- corrData[selCond,]
  DadesSel <- allData[selCond, ]  
  return(DadesSel)
}


# test
#
# (X <- round(matrix (rnorm(30)*10, ncol=6),1))
#(Y <- round(X + matrix (rnorm(30)*10, ncol=6),1))
#(rownames(X)=rownames(Y)=letters[1:nrow(X)])
#(m1<-matCorrs(X,Y))
#(m2<-matDistCorr(X,Y))
#(m12<- matAllCorrs (X, Y))
#naiveSelection(X, Y,pValCutoff=0.25)
#naiveSelection(X, Y, pValCutoff=0.25, type="Pearson")
#naiveSelection(X, Y, pValCutoff=0.25, rCutoff=0.1, type="Spearman", sortByCorrs=TRUE)
#sort1(m12,1)
#sort1(m12,3)
#dcor(X,Y,1)
#coeffRV(X,Y)
#multivCorr(X,Y)
## Adding an NA to X will make the preceeding fail
#  X[1,1] <-NA
# (m1<-matCorrs(X,Y))
# (m2<-matDistCorr(X,Y)) # yields an Error
# (m12<- matAllCorrs (X, Y)) # yields an Error
# sort1(m12,4, DEC=FALSE)
