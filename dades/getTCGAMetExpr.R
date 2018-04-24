# OBTENCIO DE DADES DEL cBioportal i LLEGIR-LES DIRECTAMENT

# Hem descarregat les dades d'expressió i metilació de cBioportal i les 
# - llegim, 
# - netegem (?) i 
# - fem concordar

x<- data_expression_merged_median_Zscores <- read.delim("gdac.20160128.cBioportal/coadread_tcga_pub/data_expression_merged_median_Zscores.txt", 
                                                        sep="\t", head = TRUE)
dim(x); class(x); head(x[1:5,1:5]); tail(x[,1:5]);
rownames(x) <- expresEntrez<- x$Entrez_Gene_Id
if(ncol(x)==225) x<- x[,-1]
# Aquesta matriu té gens i miRNAs, que cal treure.
miRNAs <- grep("hsa-*", rownames(x))
x <- x[-miRNAs,]
dim(x); class(x); head(x[1:5,1:5]); tail(x[,1:5]);

# Un cop nomes tenim Entrezs els convertim a Symbols per als gens 
require(Homo.sapiens)
anotacs<- select(Homo.sapiens, keys=rownames(x), columns= "SYMBOL", keytype="ENTREZID")
dim(anotacs); head(anotacs); tail(anotacs)
# Comptem quants gens hi ha sense símbols
sum(is.na(anotacs[,2])) 
#Treiem aquests gens i les anotacions corresponents
x<- x[!is.na(anotacs[,2]),]
anotacs<- anotacs[!is.na(anotacs[,2]),]
# Convertim els noms de les files a Symbols
rownames(x) <-anotacs[,2] 
dim(x); head(x[,1:5]); tail(x[,220:224]);

# Ens queda una matriu d'expressió amb gene Symbols a les files i amb 1094 genes que tenen al menys un NA
# Observant-los veiem que alguns tenen uns quants NAs i bastants els tene tots. Treiem aquests

quantsNAs_x <- countNAs(x)
genes2Remove <- names(quantsNAs [quantsNAs_x > 200])

expres2Keep <- x[setdiff(rownames(x),genes2Remove),]
dim(expres2Keep)
save(expres2Keep, file="TCGA-COAD-expressions0.Rda")

###########################################################################

y <- data_methylation_hm27 <- read.delim("gdac.20160128.cBioportal/coadread_tcga_pub/data_methylation_hm27.txt", 
                                        sep="\t", head = TRUE)
dim(y); head(y[,1:5]);tail(y[,1:5])

# Ens assegurem de que hi hagi un sol gen amb cada symbol

y <- y[!duplicated(y$Entrez_Gene_Id),]
y <- y[!is.na(y$Entrez_Gene_Id),]
dim(y)
# Al no estar duplicats els symbols tampoc ho estan els Entrez
y <- y[!duplicated(y$Hugo_Symbol),]
y <- y[!is.na(y$Hugo_Symbol),]
dim(y)

# Podrien els simbols ser sinonims diferents dels associats als ENtrezs?
# Per saber si es així converteixo els Entrez en Simbols i els comparo
anotacs<- select(Homo.sapiens, keys=y$Entrez_Gene_Id, columns= "SYMBOL", keytype="ENTREZID")
dim(anotacs); head(anotacs); tail(anotacs)

length(intersect(y$Hugo_Symbol,anotacs[,2]))

sum(is.na(anotacs[,2])) 
x<- x[!is.na(anotacs[,2]),]
anotacs<- anotacs[!is.na(anotacs[,2]),]

rownames(y) <- y$Hugo_Symbol
y<- y[,-c(1:2)]
head(y[1:5,1:5])

# Ens assegurem de que no quedin gens amb NAs (o amb massa NAs)

sumaY <- apply(y, 1,sum)
sum(is.na(sumaY))
yAmbNAs <- y[is.na(sumaY),]
quantsNAs <- apply(yAmbNAs, 1, function(unVec){sum(is.na(unVec))})
table(quantsNAs)
mets2Remove <- names(quantsNAs [quantsNAs > 200])

mets2Keep <- y[setdiff(rownames(y),mets2Remove),]
dim(mets2Keep)
save(mets2Keep, file="TCGA-COAD-methylations0.Rda")

##################################################################
# Podem imputar els gens amb NAs pero abans farem l'aparellament

commonGenes <- intersect(rownames(expres2Keep), rownames(mets2Keep))
length(commonGenes)
commonExpres <- expres2Keep[commonGenes,]
commonMets <- mets2Keep[commonGenes,]
commonSamples <- intersect(colnames(commonExpres), colnames(commonMets))
length(commonSamples)
commonExpres <- commonExpres[,commonSamples]
commonMets <- commonMets[, commonSamples]
commonExpres <- commonExpres[order(rownames(commonExpres)),]
commonMets <- commonMets[order(rownames(commonMets)),]

# L'ordre de les mostres és el mateix!
sum(colnames(commonMets)==colnames(commonExpres))
sum(rownames(commonMets)==rownames(commonExpres))

# revisem els NAs
sumaX <- apply(commonExpres, 1,sum); sum(is.na(sumaX))
xAmbNAs <- commonExpres[is.na(sumaX),]
quantsNAs <- apply(xAmbNAs, 1, function(unVec){sum(is.na(unVec))})
table(quantsNAs)

sumaY <- apply(commonMets, 1,sum); sum(is.na(sumaY))
yAmbNAs <- commonMets[is.na(sumaY),]
quantsNAs <- apply(yAmbNAs, 1, function(unVec){sum(is.na(unVec))})
table(quantsNAs)

# Com son pocs podem imputar-los
require(impute)
set.seed(123456)
imputed.commonExpressions <- impute.knn(as.matrix(commonExpres))
commonExpressions <- imputed.commonExpressions [["data"]]
# head(commonExpres[,1:5]); head(commonExpressions[,1:5])
set.seed(123456)
imputed.commonMethylations <- impute.knn(as.matrix(commonMets))
commonMethylations <- imputed.commonMethylations[["data"]]
# head(commonMets[,1:5]); head(commonMethylations[,1:5])

cat("TCGA Microarray data (font cBioportal) : ", dim(commonExpressions), "\n")
cat("TCGA Methylation data (font cBioportal): ", dim(commonMethylations), "\n")


# Comparació amb les dades de TCGA proporcionades pel Diego Arango

TCGAExprData <-  as.matrix(read.table(file="TCGAExprData.csv", header=TRUE, sep=",", dec=".", row.names=1))
TCGAMetilData <-  as.matrix(read.csv(file="TCGAMetilData.csv", header=TRUE, sep=",", dec=".", row.names=1))
cat("TCGA Microarray data : ", dim(TCGAExprData), "\n")
cat("TCGA Methylation data: ", dim(TCGAMetilData), "\n")
