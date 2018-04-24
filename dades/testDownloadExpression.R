require(MethylMix)
Download_GeneExpression.2 <- function (CancerSite, TargetDirectory, downloadData = TRUE)
{
  dir.create(TargetDirectory, showWarnings = FALSE)
  TCGA_acronym_uppercase = toupper(CancerSite)
  dataType = "stddata"
  if (CancerSite == "COAD") {
    dataFileTag ="COAD.Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3.2016012800.0.0.tar.gz"
  }
  cat("Searching MA data for:", CancerSite, "\\n")
  if (length(dataFileTag) == 1) {
    MAdirectories = get_firehoseData(downloadData, saveDir = TargetDirectory, 
                                     TCGA_acronym_uppercase = TCGA_acronym_uppercase, 
                                     dataFileTag = dataFileTag)
  }
  else {
    MAdirectories = c()
    for (i in 1:length(dataFileTag)) {
      MAdirectories = c(MAdirectories, get_firehoseData(downloadData, 
                                                        saveDir = TargetDirectory, TCGA_acronym_uppercase = TCGA_acronym_uppercase, 
                                                        dataFileTag = dataFileTag[i]))
    }
  }
  return(MAdirectories = MAdirectories)
}

Download_GeneExpression.2 (CancerSite="COAD", TargetDirectory=".", downloadData = TRUE)
