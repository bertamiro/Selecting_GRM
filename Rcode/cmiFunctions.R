MI <- function(x,y,h){
# x : vector of methylations
# y : vector of expressions  
# h : the std of the Gaussian kernel for density estimation 

if(length(x) != length(y))  stop("Different number of samples!")

M <- length(x) # samples
aux <- 0
two.h.square <- 2*h^2

for(i in 1:M){
  # kernel distance between the i.th data point and all other points j for each gene k
 tmpx <- x - x[i] 
 tmpx <- exp(-tmpx^2/(two.h.square))
 tmpy <- y - y[i]
 tmpy <- exp(-tmpy^2/(two.h.square))
 tmp <- tmpx*tmpy
 aux <- aux + log(M*sum(tmp)/(sum(tmpx)*sum(tmpy)))
 }
aux/M
}

cMI <- function(dataMeth, dataExp, t, h=0.3){
# dataMeth : input data methylation (a vector, not a matrix!)
# dataExp  : input data expression (a vector, not a matrix!)
# t : moves from 0 to 1
  n <- length(dataMeth)
  if(length(dataExp) != n)  stop("Different number of samples!")
  if(t < 0 | t > 1)  stop("t value is out of range")
  filter <- dataMeth < t
  ss <- sum(filter)
  if(ss != 0){
  x <- dataMeth[filter]
  y <- dataExp[filter]
  aux <- MI(x,y,h)*ss/n
  }
  else{
    aux <- 0
  }
  ss <- sum(!filter)
  if(ss != 0){
  x <- dataMeth[!filter]
  y <- dataExp[!filter]
  aux + MI(x,y,h)*ss/n
  }
  else{
    aux
  }
}

computeCMI <- function (methData, exprData){
  # check consistency
  # Provar un try/catch
  if((nrow(exprData)!=nrow(methData))||(ncol(exprData)!=ncol(methData)))
    stop("Expression amd methylation data must have same dimensions")
  tt <- seq(0,1,by=0.01)
  nt <- length(tt)
  ngenes <- nrow(exprData)
  nsamples <- ncol(exprData)
  cmi <- matrix(rep(0,nt*ngenes), ncol=nt)
  for (j in 1:nrow(methData)){
    # cat (j, " ")
    methVec <- methData[j,]
    exprVec <- exprData[j,]
    # df.xy <- na.omit(data.frame(methVec=methData[j,], exprVec=exprData[j,]))
    for(i in 1:nt){
      cmi[j ,i] <- cMI(methVec, exprVec, t=tt[i])
    }
  }
  return(cmi)
}

#'
#' cmiSelection: Wrapper to compute CMI and select genes that may be L-shaped
#'
#' @param methData
#' @param exprData
#' @param h
#' @param smallR
#' @param minCMI
#' @export
#' @example methData<- DAMetilData[1:3,]; exprData<- DAExprData[1:3,]; h=0.2; smallR=0.25; minCMI=0.1 
#' cmiSelection
#' 
cmiSelection <- function (methData, exprData, h=0.2, smallR=0.25, minCMI=0.1){
  
  # cMI function computes cMI values for different t values, from 0 to 1 and a step of 0.01 
  # Output is stored in a data frame called cMI_values.
  # For each gene, the optimal threshold is the t-value that results the minimum CMI.
  
  t <- seq(0,1,by=0.01)
  cMI_values <- matrix(NA,nrow = nrow(exprData),
                       ncol = length(t)+4,
                       byrow = TRUE)
  rownames(cMI_values) <- rownames(exprData)
  colnames(cMI_values) <- c(t,"cMI_min", "t_opt", "ratio", "meth_regulated")
  cMI_values <- as.data.frame(cMI_values)
  
  # for each row, cMI is computed for each t-value, 
  # where h is a tuning parameter for the kernel width and empirically set h = 0.2
  
  for (i in 1:nrow(methData)){
    dataMeth <- methData[i,]
    dataExp <- exprData[i,]
    for (j in 1:length(t)){
      cMI_values[i,j] <- cMI(dataMeth,dataExp,t[j], h=h)
    }
  }
  # columns cMI_min and t_opt store minimum cMI value computed above and t value used for it, respectively.
  for (i in 1:nrow(cMI_values)){
    cMI_values$cMI_min[i] <- min(cMI_values[i,1:length(t)])
  }
  for (i in 1:nrow(cMI_values)){
    cMI_values$t_opt[i] <- names(which.min(cMI_values[i,1:length(t)]))
  }
  # L-shaped genes (regulated by methylation) are selected according to three conditions 
  # (parameters were chosen according to a random permutation test [@Liu2012]).
  # 1. Ratio mincM I(t)/cM I(0) must be “small enough”, r<0.25.
  # 2. Minimum value of unconditioned CMI must be “big enough”, cMI(0)>0.1.
  # 3. Expression values must be higher in the left side of the plot than in the right side.
  for (i in 1:nrow(cMI_values)){
    cMI_values$ratio [i] <- cMI_values$cMI_min[i]/cMI_values[i,1]
  }
  mean_exp_left <- c()
  mean_exp_right <- c()
  for (i in 1:nrow(cMI_values)){
    dataMeth <- methData[i,]
    dataExp <- exprData[i,]
    filt <- dataMeth <= cMI_values$t_opt[i]
    right <- dataExp[!filt]
    left <- dataExp[filt]
    if (sum(left) == 0){
      mean_exp_left[i] <- 0
      mean_exp_right[i] <- mean(right)
    }
    if (sum(right) == 0){
      mean_exp_right[i] <- 0
      mean_exp_left[i] <- mean(left)
    }
    else{
      mean_exp_left[i] <- mean(left)
      mean_exp_right[i] <- mean(right)
    }
  }
  # Genes are labeled according three conditions above:
    for (i in 1:nrow(methData)){
      if (cMI_values$ratio[i] < smallR &
          cMI_values[i,1] > minCMI &
          mean_exp_left[i] > mean_exp_right[i]){
        cMI_values$meth_regulated[i] <- TRUE
      }
      else{
        cMI_values$meth_regulated[i] <- FALSE
      }
    }
  return(cMI_values[,c("cMI_min", "t_opt", "ratio", "meth_regulated")])
}
