> seqData <- DARNAseqData
> metData <-DAMetilData
> exprData <- DAExprData
> seqData <- DARNAseqData
> messageTitle("Scoring with default values", "=")
Scoring with default values 
=========================== 
> messageTitle("Microarrays-Methylation")
Microarrays-Methylation 
----------------------- 
> exprScores <-scoreGenesMat2(mets=metData, expres = exprData)
> selectedExprs<-exprScores[,"logicSc"]; table(selectedExprs)
selectedExprs
FALSE  TRUE 
 9140   537 
> messageTitle("RNASeq-Methylation")
RNASeq-Methylation 
------------------ 
> seqScores <-scoreGenesMat2(mets=metData, expres = seqData)
> selectedSeqs<-seqScores[,"logicSc"]; table(selectedSeqs)
selectedSeqs
FALSE  TRUE 
 8609  1068 
> messageTitle("Genes in common")
Genes in common 
--------------- 
> inCommon<- length(intersect(rownames(DAExprData[selectedExprs,]), rownames(DARNAseqData[selectedSeqs,])))
> cat(inCommon, "\n")