# TCGA data download with RTCGAToolbox
source("http://bioconductor.org/biocLite.R")
if (!(require(RTCGAToolbox))) biocLite("RTCGAToolbox")
require(RTCGAToolbox)

FIREHOSE_DATE = "20150402"
firehose_datasets = getFirehoseDatasets()
getFirehoseRunningDates(last=2)

# Colon adenocarcinoma [COAD]
tumor_type = "COAD"
COADdata <- getFirehoseData(dataset="COAD",
                            runDate="20160128",gistic2_Date="20160128",
                            mRNA_Array=TRUE)



############ on TCGAStat

library(TCGA2STAT)

ov_clinical_2 <- getTCGA(disease = "OV", data.type = "RNASeq2", clinical=TRUE)
ov_clinical_2[ov_clinical_2 == 'NA'] <- NA

clinical_variables_2 = colnames(ov_clinical_2)
clinical_variables_2_without_content = clinical_variables_2[sapply(ov_clinical_2, function(x)all(is.na(x)))]
clinical_variables_2_with_content = clinical_variables_2[! sapply(ov_clinical_2, function(x)all(is.na(x)))]

class(ov_clinical_2)
length(clinical_variables_2_with_content)

write.table(ov_clinical_2, "TCGA2STAT_ov_clinical.processed.txt", quote = F, sep = '\t')

# TCGA data download with TCGAbiolinks

library(TCGAbiolinks)

ov_clinical_3 <- TCGAquery_clinic("OV","clinical_patient")
ov_clinical_3_cleaned = ov_clinical_3[ , ! sapply(ov_clinical_3, function(x)all(x == '[Not Applicable]'))]
ov_clinical_3_cleaned[ov_clinical_3_cleaned == '[Not Available]'] <- NA

clinical_variables_3 = colnames(ov_clinical_3_cleaned)
clinical_variables_3_without_content = clinical_variables_3[sapply(ov_clinical_3_cleaned, function(x)all(is.na(x)))]
clinical_variables_3_with_content = clinical_variables_3[! sapply(ov_clinical_3_cleaned, function(x)all(is.na(x)))]

> dim (ov_clinical_3_cleaned)
[1] 587  49
> length(clinical_variables_3_with_content)
[1] 43

write.table(ov_clinical_3_cleaned, "E:/onedrive/Documents/tesis/TCGAbiolinks_ov_clinical.txt", quote = F, sep = '\t')