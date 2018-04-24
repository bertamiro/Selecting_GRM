installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(pkg)
}else{
  require(pkg, character.only=T)
  }
}

installifnot("annotate")
installifnot("GOstats")
installifnot("hwriter")
installifnot("Category")

require("annotate")
require("GOstats")
require("hwriter")
require("Category")
  
newHyperPar <- function ( geneList , geneUniverse,
              anotPackage,
              whichOntos = c("MF", "BP", "CC"),
              testDirections = c("over", "under"),
              pValueCutoff = 0.01,
              min.count = 3,
              addGeneNames = TRUE, 
              sourceOnto = "GO",
              fileName = "hyperGTestResults")
{
  hyperPar <- list(
              geneList = geneList ,
              geneUniverse = geneUniverse,
              anotPackage = anotPackage,
              whichOntos = whichOntos,
              testDirections = testDirections,
              pValueCutoff = pValueCutoff,
              min.count = min.count,
              addGeneNames = addGeneNames, 
              sourceOnto = sourceOnto,
              fileName = fileName
              )
  return(hyperPar)
}

newHiperGTestPar <- function (hyperPar){
  if (hyperPar$sourceOnto=="GO")
  {
    param = new("GOHyperGParams",
      geneIds= hyperPar$geneList,
      universeGeneIds = hyperPar$geneUniverse,
      annotation = hyperPar$anotPackage,
      ontology= hyperPar$whichOntos[1],
      pvalueCutoff= hyperPar$pValueCutoff,
      conditional=TRUE,
      testDirection=hyperPar$testDirections[1])
    }else{
    param = new("KEGGHyperGParams",
     geneIds= hyperPar$geneList,
     universeGeneIds = hyperPar$geneUniverse,
     annotation = hyperPar$anotPackage,
     pvalueCutoff = hyperPar$pValueCutoff,
     testDirection = hyperPar$testDirections[1])
  }
  return(param)
}

hiperGeometricAnalysis <- function (hyperPar, outFName="ResultsWithGenes", outputDir=".")
#
# Fa un test hipergeometric tant per KEGG com per GO segons els parametres que li passem
#
{
# Crea l'hiperparàmetre

  hGparam <- newHiperGTestPar(hyperPar)
  sourceOnto <- hyperPar$sourceOnto
  
  if (sourceOnto=="GO")
  {
     ontos <- hyperPar$whichOntos
  }else{
     ontos <- "KEGG"
  }
  directions <- hyperPar$testDirections

# Una funcio per completar les anotacions dels resultats
     anotPack <- annotation (hGparam)
     getSymbol <- function (x) {
       if (length(x)>0){
         simbols <- getSYMBOL(x, anotPack)
       }else{
         simbols <- NULL
       }
     return(simbols)
     }
  
# Fa el test hiperGeometric

  for (i in 1:length(ontos)){
     if (sourceOnto=="GO")
        {ontology(hGparam) <- ontos[i]}

     hyperRes <- hyperGTest(hGparam)

     # Resum dels resultats
 
     sumari <- summary(hyperRes, p=pvalueCutoff(hGparam))
     
     # Informe en HTML
     if (sourceOnto=="GO"){
       fName <- paste(hyperPar$fileName,sourceOnto, ontos[i], sep=".")
     }else{
       fName <- paste(hyperPar$fileName,sourceOnto,sep=".")
     }
     htmlReport(hyperRes, file = paste(fName, "html", sep="."), summary.args=list("htmlLinks"=TRUE))
     if ( hyperPar$addGeneNames){
        # Repeat summaries annnotating categories with the corresponding genes
       genes <-  lapply(geneIdsByCategory(hyperRes), getSymbol)
       genes2<- data.frame(genes=unlist(lapply(genes, FUN = function (x) paste (x,collapse=" "))))
       # Un truquet per concatenar-los
       sum2 <- summary(hyperRes, p=1)
       genes2 <- genes2[rownames(genes2) %in% as.character(sum2[,1]),]
       # Retoquem el sumari per poder enganxar les diverses ontologies
       Report <-cbind(onto=rep(ontos[i],nrow(sum2)),sum2, genes2)
       colnames(Report)[2]<-"ID"
     }else{
       Report <-cbind(onto=rep(ontos[i],nrow(sum2)),sum2)
       colnames(Report)[2]<-"ID"
     }
     # Cal controlar que passa quan no es selecciona res
     if(i==1){
       ReportSig <- Report[1:nrow(sumari),]
     }else{
       ReportSig <- rbind(ReportSig,  Report[1:nrow(sumari),])
     }
     
   }
  # Substituir per una crida al paquet que escriu taules en html
  
  txtFName <- paste(paste(hyperPar$fileName, sourceOnto, sep="."), outFName, "txt", sep=".")
  write.table(ReportSig, file.path(outputDir, txtFName), sep="\t")
  htmlFName <- paste(paste(hyperPar$fileName, sourceOnto, sep="."), outFName, "html", sep=".")
  hwrite(ReportSig,file.path(outputDir, htmlFName)  )

  return(ReportSig)
}
