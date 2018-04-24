# PART 2: DESCARREGAR I PROCESSAR DADES DE METILACIO I D'EXPRESIO AMB METHYLMIX
require(MethylMix)
cancerSite <- "COAD"
targetDirectory <- paste0(getwd(), "/")

# 1.2 Downloading methylation data
# METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory, TRUE)
# Els baixo manualment i elimino el Met450k
# MetDir27k <- "./gdac_20160128/gdac.broadinstitute.org_COAD.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/"
MetDir27k <- paste0(getwd(), "/","gdac_20160128/Met27kFiles")

METdirectories <- list(METdirectory27k=MetDir27k)

# Processing methylation data
METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)

# Saving methylation processed data
saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
# Clustering methylation data
res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])
# Saving methylation clustered data
toSave <- list(METcancer = res[[1]], METnormal = res[[2]], ProbeMapping = res$ProbeMapping)
saveRDS(toSave, file = paste0(targetDirectory, "MET_", cancerSite, "_Clustered.rds"))

